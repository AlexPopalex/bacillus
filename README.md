This is the repository for the scripts used in the manuscript

### Evolution of alternative reproductive systems in Bacillus stick insects

by Guillaume Lavanchy*, Alexander Brandt*, Marc Bastardot, Zoé Dumas, Marjorie Labédan, William Toubiana, Patrick Tran Van, Andrea Luchetti, Valerio Scali, Barbara Mantovani, Tanja Schwander.



# mtDNA phylogeny
We added one reference individual per haplotype from this paper: https://doi.org/10.1006/mpev.2000.0850
We aligned the sequences with muscle 3.8.1551:
```
muscle -in All_filtered_with_refs.fa -phys > phylip.tmp
NSEQ=$( head -1 phylip.tmp | cut -f 1 -d ' ' )
NSITES=$( head -1 phylip.tmp | cut -f 2 -d ' ' )
echo $NSEQ $NSITES > All_filtered_with_refs_aligned.phylip
tail -n +2 phylip.tmp | tr -d '\n' | tr '@' '\n' | sed 's/ /@/1' | tr -d ' ' | tr '@' ' ' >> All_filtered_with_refs_aligned.phylip
```

The resulting alignment is available in this repository: `All_filtered_with_refs_aligned.phylip`. We then reconstructed the phylogeny with iqtree 2.0.6:
```
iqtree2 -T 1 -s Sicily_with_refs_aligned.phylip -m MFP -B 1000
```
The resulting tree was rerooted manually in FigTree using the outgroup (*Clonopsis*) and analysed in R along with the rest. The rerooted version is available in this repository: `All_filtered_with_refs_aligned.phylip_rerooted.newick`.

# Population structure
We assigned each individual to a lineage and assessed the genomic composition of individuals mostly in `R` v. 4.1.1. The code used for the analyses is given in **Bacillus_analyses.R**.

We ran `admixture` on a subset of autosomal loci. We filtered the original vcf file with `vcftools`:
```bash
vcftools --vcf Brsri_v3_gatk.vcf \
--recode-INFO-all \
--not-chr Brsri_v3_scf2 \
--mac 3 \
--max-alleles 2 \
--maxDP 40 \
--min-alleles 2 \
--minDP 8 \
--minQ 20 \
--max-missing 0.5 \
--out Sicily_filtered_shared_gatk_an23_ql20_snps_dp08_mc3_miss50 \
--positions positions_genotyped_in_all_sspp.txt \
--recode \
--remove-indels
```

We then ran admixture over 10 replicates using the following script:
```bash
### This script takes 3 arguments:
### $1 is the input file in vcf format. The chromosome names must be numeric values because plink is dumb.
### $2 is the number of K that we want to test.
### $3 is the number of replicates that we want to run.
### Usage example: sbatch admixture.sh cleaned.vcf 6 10

# pretending all markers are on chromosome 1
grep '#' $1 > headers.tmp
grep -v '#' $1 | cut -f 2- > vcf2.tmp
NLINES=$(wc -l vcf2.tmp | cut -f 1 -d ' ')
yes 1 | head -n $NLINES > vcf1.tmp
paste vcf1.tmp vcf2.tmp > vcf3.tmp
cat headers.tmp vcf3.tmp > cleaned.vcf
rm *.tmp

module load gcc/9.3.0
module load plink-ng/2.0.20200727

plink2 --vcf cleaned.vcf --make-bed --out cleaned --allow-extra-chr --chr-set 28

for i in $(seq 1 $3)
do
mkdir rep${i}
cd rep${i}

~/software/admixture -s time --supervised --cv=10 ../cleaned.bed $2 &> log_K$2.out

grep -h CV log*.out > CV_results_rep${i}.txt

cd ../

done

# Summarizing the CV error

cut -f 1 -d ')' rep1/CV_results_rep1.txt | cut -f 2 -d '=' > k.tmp
for j in $(ls rep*/CV_results_rep*.txt | cut -f 1 -d '/')
do
cat ${j}/CV_results_${j}.txt | cut -f 4 -d ' ' > ./${j}.tmp
done
paste k.tmp rep*.tmp > CV_summary_across_rep.txt
rm *.tmp
```


