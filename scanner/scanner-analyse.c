/* See COPYING file for copyright and license details. */

#include "scanner-dump.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "parse_args.h"
#include "scanner-common.h"
#include "nproc.h"

extern gboolean verbose;
static gboolean pkLevels = FALSE;
static double interval = 1.0;
static double inputscale = 1.0;
extern gchar *decode_to_file;

static GOptionEntry entries[] =
{
    { "interval", 'i', 0, G_OPTION_ARG_DOUBLE, &interval, NULL, NULL },
    { "scale",    's', 1, G_OPTION_ARG_DOUBLE, &inputscale, NULL, NULL },
    { "pklevels", 'l', 0, G_OPTION_ARG_NONE, &pkLevels, NULL, NULL },
    { NULL, 0, 0, G_OPTION_ARG_NONE, NULL, NULL, 0 }
};

static ebur128_state *st;
static double samplePeakLast[8], samplePeakAll[8], truePeakLast[8], truePeakAll[8];
static double samplePeakLastMax[8], truePeakLastMax[8];
static double loudness_m, loudness_s, loudness_i;
static double loudness_m_max, loudness_s_max;
static double loudness_lu, rel_thr;


static double level( double v )
{
    if ( v < 1E-20 )
        return -200.0;
    else
        return 20.0 * log10( v );
}

static void resetLevelMax(const unsigned nchans_max)
{
    for (unsigned ch = 0; ch < nchans_max; ++ch)
        samplePeakLastMax[ch] = truePeakLastMax[ch] = 0.0;
    loudness_m_max = loudness_s_max = -DBL_MAX;
}

static void updateLevels(const unsigned nchans_max)
{
    ebur128_loudness_momentary(st, &loudness_m);
    ebur128_loudness_shortterm(st, &loudness_s);
    ebur128_loudness_global(st, &loudness_i);
    ebur128_loudness_range(st, &loudness_lu);
    ebur128_relative_threshold(st, &rel_thr);

    for (unsigned ch = 0; ch < nchans_max; ++ch) {
        ebur128_prev_sample_peak(st, ch, &samplePeakLast[ch]);
        ebur128_sample_peak(st, ch, &samplePeakAll[ch]);
        ebur128_prev_true_peak(st, ch, &truePeakLast[ch]);
        ebur128_true_peak(st, ch, &truePeakAll[ch]);
        
        if ( samplePeakLastMax[ch] < samplePeakLast[ch] )
            samplePeakLastMax[ch] = samplePeakLast[ch];
            
        if ( truePeakLastMax[ch] < truePeakLast[ch] )
            truePeakLastMax[ch] = truePeakLast[ch];        
    }

    if ( loudness_m_max < loudness_m )
        loudness_m_max = loudness_m;
    if ( loudness_s_max < loudness_s )
        loudness_s_max = loudness_s;
}


static void dump_loudness_info(struct filename_list_node *fln, int *ret)
{
    struct input_ops* ops = NULL;
    struct input_handle* ih = NULL;
    float *buffer = NULL;

    int result;
    static size_t nr_frames_read;
    static size_t frames_counter, frames_needed;
    const int r128_mode = EBUR128_MODE_M | EBUR128_MODE_S | EBUR128_MODE_I
        | EBUR128_MODE_LRA | EBUR128_MODE_SAMPLE_PEAK | EBUR128_MODE_TRUE_PEAK;

    result = open_plugin(fln->fr->raw, fln->fr->display, &ops, &ih);
    if (result) {
        *ret = EXIT_FAILURE;
        goto free;
    }

    const unsigned nchans = ops->get_channels(ih);
    const long srate = ops->get_samplerate(ih);
    if (!st) {
        st = ebur128_init(nchans, srate, r128_mode);
        if (!st) abort();
    } else {
        if (!ebur128_change_parameters(st, nchans, srate)) {
            frames_counter = 0;
        }
    }

    result = ops->allocate_buffer(ih);
    if (result) abort();
    buffer = ops->get_buffer(ih);

    frames_needed = (size_t) (interval * (double) st->samplerate + 0.5);
    const unsigned nchans_max = ( nchans > 8 ) ? 8 : nchans;

    unsigned long frame_i_counter = 0;

    printf("time/sec"
           "\tmomentary(400ms)/LUFS"
           "\tshortterm(3sec)/LUFS"
           "\tglobal_integrated/LUFS"
           "\t"
           "\trange_LRA/LU"
           "\trelative_threshold/LUFS"
           );
    for (unsigned ch = 0; ch < nchans_max; ++ch)
        printf("\t\tlast_SPK_%d\tall_SPK_%d"
               "\tlast_TPK_%d\tall_TPK_%d", ch, ch, ch, ch );
    printf("\n");


    while ((nr_frames_read = ops->read_frames(ih))) {
        float* tmp_buffer = buffer;
        while (nr_frames_read > 0) {
            if (frames_counter + nr_frames_read >= frames_needed) {
                const double t = (double)frame_i_counter / (double)srate;
                frame_i_counter += frames_counter + nr_frames_read;
                if (inputscale != 1.0) {
                    float scale = inputscale;
                    for (size_t k = 0; k < (frames_needed - frames_counter) * nchans_max; ++k )
                        tmp_buffer[k] *= scale;
                }
                result = ebur128_add_frames_float(st, tmp_buffer,
                                                  frames_needed - frames_counter);
                if (result) abort();
                tmp_buffer += (frames_needed - frames_counter) * st->channels;
                nr_frames_read -= frames_needed - frames_counter;
                frames_counter = 0;

                updateLevels(nchans_max);

                printf("%.1f"
                       "\t%.1f\t%.1f\t%.1f"
                       "\t\t%.1f\t%.1f", t, loudness_m_max, loudness_s_max, loudness_i, loudness_lu, rel_thr);
                if ( pkLevels )
                {
                    for (unsigned ch = 0; ch < nchans_max; ++ch)
                        printf("\t\t%.1f\t%.1f\t%.1f\t%.1f"
                            , level(samplePeakLastMax[ch]), level(samplePeakAll[ch])
                            , level(truePeakLastMax[ch]), level(truePeakAll[ch]) );
                }
                else
                {
                    for (unsigned ch = 0; ch < nchans_max; ++ch)
                        printf("\t\t%6.4f\t%6.4f\t%6.4f\t%6.4f"
                            , samplePeakLastMax[ch], samplePeakAll[ch]
                            , truePeakLastMax[ch], truePeakAll[ch] );
                }
                printf("\n");
                
                resetLevelMax(nchans_max);

                
            } else {
                if (inputscale != 1.0) {
                    float scale = inputscale;
                    for (size_t k = 0; k < nr_frames_read * nchans_max; ++k )
                        tmp_buffer[k] *= scale;
                }
                result = ebur128_add_frames_float(st, tmp_buffer, nr_frames_read);
                if (result) abort();
                updateLevels(nchans_max);
                tmp_buffer += (nr_frames_read) * st->channels;
                frames_counter += nr_frames_read;
                nr_frames_read = 0;
            }
        }
    }

  free:
    if (ih) ops->free_buffer(ih);
    if (!result) ops->close_file(ih);
    if (ih) ops->handle_destroy(&ih);
}

int loudness_analyse(GSList *files)
{
    int ret = 0;

    g_slist_foreach(files, (GFunc) dump_loudness_info, &ret);
    if (st) ebur128_destroy(&st);

    return ret;
}

gboolean loudness_analyse_parse(int *argc, char **argv[])
{
    if (decode_to_file) {
        fprintf(stderr, "Cannot decode to file in dump mode\n");
        return FALSE;
    }

    if (!parse_mode_args(argc, argv, entries)) {
        if (*argc == 1) fprintf(stderr, "Missing arguments\n");
        return FALSE;
    }

    return TRUE;
}
