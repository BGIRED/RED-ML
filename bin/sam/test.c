#include "stdio.h"
#include "bam.h"
#include "bgzf.h"
#include "faidx.h"
#include "knetfile.h"
#include "glf.h"
#include "sam.h"
#include "kstring.h"
#include "razf.h"
#include "sam_header.h"


int main()
{	
	samfile_t *in = 0, *out = 0;
        int slx2sngr = 0;
        char in_mode[5], out_mode[5], *fn_out = 0, *fn_list = 0, *fn_ref = 0;
        strcpy(in_mode, "r"); strcpy(out_mode, "w");
        strcat(in_mode, "b");
                printf("start******************************************************\n");
        char *filename = "/ifs1/RD/shaohaojing/DAI/1000306/1000306_HUMgqbRLJDIAAPE_091015_I58_FC42UC6AAXX_L2_HUMgqbRLJDIAAPE_1.fq.bam.sort.bam";
        if ((in = samopen(filename, in_mode, fn_list)) == 0)
        {
                fprintf(stderr, "[main_samview] fail to open file for reading.\n");
                exit(0);
        }
        if (in->header == 0) {
                fprintf(stderr, "[main_samview] fail to read the header.\n");
                exit(0);
        }
        if ((out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0) {
                fprintf(stderr, "[main_samview] fail to open file for writing.\n");
                exit(0);
        }

        bam1_t *b = bam_init1();
        int r;
        char *s;
        while((r = samread(in, b)) >= 0)
        {
                if (!__g_skip_aln(in->header, b)) {
                        s = bam_format1_core(out->header, b, out->type>>2&3);
                }
                printf("%s\n", s);
                free(s);
//                printf("%s", b->data);
//                printf("\n++++++++++%d+++++++++++", r);
//              printf("********************\n");
        }
        bam_destroy1(b);
        samclose(in);
        samclose(out);
}

