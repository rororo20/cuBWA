
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "smems/bntseq.h"
#include "smems/bwa.h"
#include "smems/bwt.h"
#include "smems/FMI_search.h"
#include "smems/macro.h"
#include "string.h"
#include "smems/utils.h"

uint64_t proc_freq, tprof[LIM_R][LIM_C], prof[LIM_R];

int bwa_index(int argc, char *argv[])  // the "index" command
{
    int c;
    char *prefix = 0;
    bool order = false;
    while ((c = getopt(argc, argv, "p:r")) >= 0) {
        if (c == 'p')
            prefix = optarg;
        else if (c == 'r')
            order = true;
        else
            return 1;
    }

    if (optind + 1 > argc) {
        fprintf(stderr, "Usage: bwa-mem2 index [-p prefix] [-r] <in.fasta>\n");
        return 1;
    }
    if (prefix == 0) prefix = argv[optind];
    bwa_idx_build(argv[optind], prefix, order);
    return 0;
}

int bwa_idx_build(const char *fa, const char *prefix, bool order)
{
    extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

    clock_t t;
    int64_t l_pac;

    {  // nucleotide indexing
        gzFile fp = xzopen(fa, "r");
        t = clock();
        fprintf(stderr, "[bwa_index] Pack FASTA... ");
        l_pac = bns_fasta2bntseq(fp, prefix, 1);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        err_gzclose(fp);
        FMI_search *fmi = new FMI_search(prefix, order);
        fmi->build_index();
        delete fmi;
    }
    return 0;
}
int main(int argc, char *argv[])
{
    uint64_t tim = __rdtsc();
    sleep(1);
    proc_freq = __rdtsc() - tim;

    int ret = -1;
    if (argc < 2) {
        fprintf(stderr, "cuBWA tools:\n");
        return -1;
    }
    if (strcmp(argv[1], "index") == 0) {
        uint64_t tim = __rdtsc();
        ret = bwa_index(argc - 1, argv + 1);
        fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim) * 1.0 / proc_freq);
        return ret;
    }
    else if (strcmp(argv[1], "mem") == 0) {
        return ret;
    }
    return 0;
}