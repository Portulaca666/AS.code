###############################run in linux########################################
####single-end RNA-seq analysis
pkurun-fat4way 1 10 prefetch SRR7107{499..517}
##转换为fastq
for id in SRR7107{499..517}
do fastq-dump ${id}/${id}.sra --split-files -O /home/xjhuang_pkuhpc/gpfs1/help/LQQ/raw.fastq/
  done

fastq-dump SRR7107518/SRR7107518.sra --split-files -O /home/xjhuang_pkuhpc/gpfs1/help/LQQ/raw.fastq/
  
  out_dir=/home/xjhuang_pkuhpc/gpfs1/help/LQQ/clean
for id in SRR7107{499..517}
do trim_galore  --quality 20 -a AGATCGGAAGAGC --length 20 --cores 9 -o ${out_dir}  /home/xjhuang_pkuhpc/gpfs1/help/LQQ/raw.fastq/${id}_1.fastq
done

out_dir=/home/xjhuang_pkuhpc/gpfs1/help/LQQ/fastqc
fastqc -o out_dir  -t 10 input.fq

##run   rMATs
s1_path=/home/xjhuang_pkuhpc/gpfs1/help/LQQ/clean/trimed_fastq/ss1.fastq.csv
s2_path=/home/xjhuang_pkuhpc/gpfs1/help/LQQ/clean/trimed_fastq/ss2.fastq.csv

source /appsnew/source/rmats_turbo_v4.2.0.sh
python /appsnew/usr/python/Miniconda/miniconda3-py311_23.5.2/envs/rmats_turbo_v4_2_0/rmats.py -t single --s1 ${s1_path} --s2 ${s2_path}  --gtf /home/xjhuang_pkuhpc/lustre2/genome/NCBI/hg38/hg38.ncbiRefSeq.gtf --bi /home/xjhuang_pkuhpc/lustre2/genome/STAR/hg38  --readLength 75 --variable-read-length    --nthread 18 --od /home/xjhuang_pkuhpc/gpfs1/help/LQQ/rmats/output --tmp /home/xjhuang_pkuhpc/gpfs1/help/LQQ/rmats/tmp_output


###stats events number
for i in RI SE A5SS A3SS MXE; do 
awk -F "\t" '{print $0}' ${i}.MATS.JCEC.txt |wc -l; 
done
##stats genes number
for i in RI SE A5SS A3SS MXE; do 
awk -F "\t" '{print $0}' ${i}.MATS.JCEC.txt |cut -f2|sort|uniq |wc -l;
done

##run   rmats2sashimiplot 
sed -n '/geneSymbol\|Fn1/p'   SE.MATS.JC.txt >  grep_Fn1_SE.MATS.JC.txt
sed -n '/geneSymbol\|Plod2/p'   SE.MATS.JC.txt >  grep_Plod2_SE.MATS.JC.txt
test_coordinate_output=/home/xjhuang_pkuhpc/gpfs1/help/LQQ/mouse/RNA_mouse/rMATs/NC_KO_C/rmats2sashimiplot_output
output_paths=/home/xjhuang_pkuhpc/gpfs1/help/LQQ/mouse/RNA_mouse/rMATs/NC_KO_C/output
rmats2sashimiplot --b1 `cat /home/xjhuang_pkuhpc/gpfs1/help/LQQ/mouse/RNA_mouse/rMATs/NC_KO/tmp_output/sample_WT` --b2 `cat /home/xjhuang_pkuhpc/gpfs1/help/LQQ/mouse/RNA_mouse/rMATs/NC_KO/tmp_output/sample_KO`  --l1 SampleWT --l2 SampleKO --event-type SE  -e ${output_paths}/grep_Plod2_SE.MATS.JC.txt  -o ${test_coordinate_output}
rmats2sashimiplot --b1 `cat /home/xjhuang_pkuhpc/gpfs1/help/LQQ/mouse/RNA_mouse/rMATs/NC_KO/tmp_output/sample_WT` --b2 `cat /home/xjhuang_pkuhpc/gpfs1/help/LQQ/mouse/RNA_mouse/rMATs/NC_KO/tmp_output/sample_KO`  --l1 SampleWT --l2 SampleKO --event-type SE  -e ${output_paths}/grep_Fn1_SE.MATS.JC.txt  -o ${test_coordinate_output}
