/* See COPYING file for copyright and license details. */

#ifndef SCANNER_ANALYSE_H
#define SCANNER_ANALYSE_H

#include <glib.h>

int loudness_analyse(GSList *files);
gboolean loudness_analyse_parse(int *argc, char **argv[]);

#endif /* end of include guard: SCANNER_ANALYSE_H */
