###############################################################################
## Pipeline of PiRNA cluster extraction from GRCh37 genome version of 1000GP ##
###############################################################################
(0) Preliminary steps:

If you do not have administrator privileges, do the following:

mkdir ~/local
mkdir ~/local

Go to your .bashrc profile and add the following lines:

nano ~/.bashrc
PATH=~/Scripts:$PATH # Only to run our personalized scripts
PATH=~/local/bin:$PATH
export X # We can add any environment variables required, such as BCFTOOLS_PLUGINS.

(1) Installation of key features:

|A| Bcftools: This pipeline uses the experimental version because it has a feature to replace ./. with 0/0 in merge.

# Install the regular version [Not necessary]

wget https://github.com/samtools/htslib/releases/download/1.3/bcftools-1.3.tar.bz2
tar -xvf bcftools-1.3.tar.bz2
cd bcftools-1.3
make
make prefix=~/local install

# Install the experimental version:

git clone --branch=develop --recursive git://github.com/pd3/bcftools.git
cd bcftools; make

##############################################################
## PIRNA CLUSTER EXTRACTIONS FROM GENOTYPES OF 1KGP, GRCh37 ##
##############################################################

#(1) Extract all fragments from each chromosome and index them with TABIX, storing in the corresponding to each chromosome directory

CHRS=`seq 1 22`

for CHR in $CHRS; do
  mkdir -p chr${CHR}
  POSITIONS=`cat cluster_positions_GRCh37/Human_clusters_hg19_CHR${CHR}.bed | sed "s/^chr\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1:\2-\3/g" | tr "\n" " "`
  tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${POSITIONS} >chr${CHR}/chr${CHR}.vcf
  bcftools view -O -S ^InbredIndividuals.txt chr${CHR}/chr${CHR}.vcf > chr${CHR}/chr${CHR}_filtered.vcf
  vcf-sort chr${CHR}/chr${CHR}_filtered.vcf > chr${CHR}/chr${CHR}_sort.vcf
  bgzip chr${CHR}/chr${CHR}_sort.vcf
  tabix -p vcf chr${CHR}/chr${CHR}_sort.vcf.gz
done


#########################################################
## PIRNA CLUSTER EXTRACTIONS FROM CHIMPANZEE ALIGNMENT ##
#########################################################

#(1) Download MFA files from Vista Browser, fasta reference files from UCSC browser, convert human-chimp alingment into vcf file with MFAtoVCF.py script

mkdir -p Outgroup_GRCh37

CHRS=`seq 1 22`

for CHR in $CHRS; do
    mkdir -p chr${CHR}
    wget http://pipeline.lbl.gov/data/hg19_panTro4/chr${CHR}.mfa.gz > chr${CHR}/chr${CHR}.mfa.gz
	gunzip chr${CHR}/chr${CHR}.mfa.gz
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${CHR}.fa.gz > chr${CHR}/chr${CHR}.fa.gz
    gunzip chr${CHR}/chr${CHR}.fa.gz
	grep -v "score" chr${CHR}/chr${CHR}.mfa > chr${CHR}/chr${CHR}.wscore.mfa
    /home/path/to/soft/Python-3.5.0/python MFAtoVCF.py -s chr${CHR}/chr${CHR}.fa -q Chimp -c ${CHR} chr${CHR}/chr${CHR}.wscore.mfa 
done

#the same for chromosome X
mkdir -p chrX
    wget http://pipeline.lbl.gov/data/hg19_panTro4/chrX.mfa.gz > chrX/chrX.mfa.gz
	gunzip chrX/chrX.mfa.gz
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz > chrX/chrX.fa.gz
    gunzip chrX/chrX.fa.gz
	grep -v "score" chrX/chrX.mfa > chrX/chrX.wscore.mfa
    /home/path/to/soft/Python-3.5.0/python chrX/MFAtoVCF.py -s chrX/chrX.fa -q Chimp -c X chrX/chrX.wscore.mfa 
	
#(2) Extract all fragments from each chromosome and index them with TABIX, storing in the corresponding to each chromosome directory

CHRS=`seq 1 22`

for CHR in $CHRS; do
  mkdir -p chr${CHR}
  POSITIONS=`cat cluster_positions_GRCh37/Human_clusters_hg19_CHR${CHR}.bed | sed "s/^chr\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1:\2-\3/g" | tr "\n" " "`
  tabix -f -h chr${CHR}/chr${CHR}.vcf.gz ${POSITIONS} > chr${CHR}/chr${CHR}_chimp.vcf

  vcf-sort chr${CHR}/chr${CHR}_chimp.vcf > chr${CHR}/chr${CHR}_chimp_sort.vcf
  bgzip chr${CHR}/chr${CHR}_chimp_sort.vcf
  tabix -p vcf chr${CHR}/chr${CHR}_chimp_sort.vcf.gz
done

#Test on chromosome X 

  POSITIONS=`cat cluster_positions_GRCh37/Human_clusters_hg19_CHRX.bed | sed "s/^chr\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1:\2-\3/g" | tr "\n" " "`
  tabix -f -h chrX/chrX.vcf.gz ${POSITIONS} > chrX/chrX_chimp.vcf
  vcf-sort chrX/chrX_chimp.vcf > chrX/chrX_chimp_sort.vcf
  bgzip chrX/chrX_chimp_sort.vcf
  tabix -p vcf chrX/chrX_chimp_sort.vcf.gz

#######################################################
## MERGE VCF.GZ FILES OF HUMAN GENOMES WITH OUTGROUP ##
#######################################################

mkdir -p Merged_chromosomes
CHRS=`seq 1 22`

for CHR in $CHRS; do
$bcf/bcftools merge  -Oz --missing-to-ref path/to/human/chr${CHR}_sort.vcf.gz path/to/chimp/chr${CHR}_chimp_sort.vcf.gz > path/to/Merged_chromosomes/chr${CHR}_merged.vcf.gz
tabix -p vcf path/to/Merged_chromosomes/chr${CHR}_merged.vcf.gz;
done

#Merge chrX for test
$bcf/bcftools merge  -Oz --missing-to-ref path/to/human/chrX_sort.vcf.gz path/to/chimp/chrX_chimp.vcf.gz > path/to/Merged_chromosomes/chrX_merged.vcf.gz
tabix -p vcf path/to/Merged_chromosomes/chrX_merged.vcf.gz


###############################################
## EXTRACTION OF INTERGENIC REGIONS FOR  MKT ##
###############################################

#(1) Obtaining the coordinates of intergenic regions

#Download the GENCODE annotations for GRCh37/hg19 genome reference version

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz

#To define intronic regions, we need to define the gene i.e. obtain the gene coordinates

zcat gencode.v25lift37.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5}' > gencode_GRCh37_gene.bed

/home/path/to/software/bedtools/bin/sortBed -i gencode_GRCh37_gene.bed > gencode_GRCh37_gene_temp.bed

mv -f gencode_GRCh37_gene_temp.bed gencode_GRCh37_gene.bed

#And finally to define intergenic regions, we use complementBed to find regions not covered by genes. 
#To create a hg19_chrom_info.txt file, use the fetchChromSizes executable available at http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
#to create the hg19_chrom_info.txt file for hg19.

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes 

chmod +x fetchChromSizes

./fetchChromSizes hg19 > hg19_chrom_info.txt

/home/rpath/to/software/bedtools/bin/complementBed -i gencode_GRCh37_gene.bed -g hg19_chrom_info.txt > gencode_GRCh37_intergenic.bed

###################################################################
##Defining overlapping regions with genes, pseudogenes and polyAs##
###################################################################

/home/path/to/software/bedtools/bin/intersectBed -a Clusters_GRCh37.bed -b gencode_GRCh37_gene.bed -wo > clusters_overlap_genes.bed

zcat gencode.v25.2wayconspseudos.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="transcript" {print $1,$4,$5}' > gencode_GRCh37_pseudogenes.bed

/home/path/to/software/bedtools/bin/intersectBed -a Clusters_GRCh37.bed -b gencode_GRCh37_pseudogenes.bed -wo > clusters_overlap_pseudogenes.bed

olga@olga-HP-Compaq-8200-Elite-MT-PC:~/Downloads$ zcat gencode.v25.polyAs.gtf.gz | awk 'BEGIN{OFS="\t";} {print $1,$4,$5}' > gencode_GRCh37_polyAs.bed

odolgova@andromeda:~/Files$ /home/path/to/software/bedtools/bin/intersectBed -a Clusters_GRCh37.bed -b gencode_GRCh37_polyAs.bed -wo > clusters_overlap_polyAs.bed

#Concatenate the intergenic regions in BED file#
odolgova@andromeda:~/Files/Closest_intergenic_regions/Human_intergenic/Intergenic_regions_hg19$ cat Intergenic_hg19_CHR*.bed | sort -k 1,1 -k2,2n > Intergenic_hg19_all.bed

###############Definining overlapped with intergenic regions########################

/home/raquel/pj15011/software/bedtools/bin/intersectBed -a Intergenic_hg19_all.bed -b gencode_GRCh37_pseudogenes.bed -wo > intergenic_overlap_pseudogenes.bed

/home/raquel/pj15011/software/bedtools/bin/intersectBed -a Intergenic_hg19_all.bed -b gencode_GRCh37_polyAs.bed -wo > intergenic_overlap_polyAs.bed

zcat gencode.v25.tRNAs.gtf.gz | awk 'BEGIN{OFS="\t";} {print $1,$4,$5}' > gencode_GRCh37_tRNAs.bed

/home/raquel/pj15011/software/bedtools/bin/intersectBed -a Intergenic_hg19_all.bed -b gencode_GRCh37_tRNAs.bed -wo > intergenic_overlap_tRNAs.bed


#######################################################################################
##Discarding the regions overlapping with genes, pseudogenes and polyAs from BED files#
#######################################################################################

odolgova@andromeda:~/Files/Merged_GRCh37/overlapping$ cat intergenic_overlap_*.bed | sort  -k 1,1 -k2,2n > intergenic_overlap_all.bed
odolgova@andromeda:~/Files/Merged_GRCh37/overlapping$ cat clusters_overlap_*.bed | sort  -k 1,1 -k2,2n > clusters_overlap_all.bed

##############################################################
##Substraction of the gene/pseudogene regions from BED files##
##############################################################

/home/path/to/software/bedtools/bin/subtractBed -a Clusters_GRCh37.bed -b clusters_overlap.bed SV_map.bed > new_clusters.bed

/home/path/to/software/bedtools/bin/subtractBed -a Intergenic_hg19_all.bed -b intergenic_overlap.bed SV_map.bed > new_intergenics.bed


#############################################################
##Extraction of superpopulations from chromosomic vcf files##
#############################################################

#For clusters
CHRS=`seq 1 22`
SETS=`seq 1 5`
for CHR in $CHRS; do
echo "starting chr${CHR}"
 rm chr${CHR}.vcf.gz.tbi
 gunzip chr${CHR}.vcf.gz
for SET in $SETS; do
echo "starting population ${SET}"
 vcf-subset -c ${SET}_individuals.txt chr${CHR}.vcf > chr${CHR}_${SET}.vcf
bgzip chr${CHR}_${SET}.vcf
tabix -p vcf chr${CHR}_${SET}.vcf.gz
done
bgzip chr${CHR}.vcf
tabix -p vcf chr${CHR}.vcf.gz
done

#test on X chromosome
SETS=`seq 1 5`

rm chrX.vcf.gz.tbi
gunzip chrX.vcf.gz
for SET in $SETS; do
vcf-subset -c ${SET}_individuals.txt chrX.vcf > chrX_${SET}.vcf
bgzip chrX_${SET}.vcf
tabix -p vcf chrX_${SET}.vcf.gz
done
bgzip chrX.vcf
tabix -p vcf chrX.vcf.gz

vcf-subset -c 5_individuals.txt chrX.vcf > chrX_5.vcf
bgzip chrX_5.vcf
bgzip chrX.vcf
tabix -p vcf chrX_5.vcf.gz
tabix -p vcf chrX.vcf.gz

#For intergenic regions
CHRS=`seq 1 22`
SETS=`seq 1 5`
for CHR in $CHRS; do
echo "starting chr${CHR}"
rm chr${CHR}_intergenic.vcf.gz.tbi
gunzip chr${CHR}_intergenic.vcf.gz
for SET in $SETS; do
echo "starting population ${SET}"
vcf-subset -c ${SET}_individuals.txt chr${CHR}_intergenic.vcf > chr${CHR}_${SET}_intergenic.vcf
bgzip chr${CHR}_${SET}_intergenic.vcf
tabix -p vcf chr${CHR}_${SET}_intergenic.vcf.gz
done
bgzip chr${CHR}_intergenic.vcf
tabix -p vcf chr${CHR}_intergenic.vcf.gz
done

#test on X chromosome
SETS=`seq 1 5`

rm chrX_intergenic.vcf.gz.tbi
gunzip chrX_intergenic.vcf.gz
for SET in $SETS; do
vcf-subset -c ${SET}_individuals.txt chrX_intergenic.vcf > chrX_${SET}_intergenic.vcf
bgzip chrX_${SET}_intergenic.vcf
tabix -p vcf chrX_${SET}_intergenic.vcf.gz
done
bgzip chrX_intergenic.vcf
tabix -p vcf chrX_intergenic.vcf.gz


###################################################################
##Substraction of the gene/pseudogene regions for each population##
###################################################################

/home/path/to/software/bedtools/bin/subtractBed -a Clusters_GRCh37.bed -b clusters_overlap.bed SV_map.bed > new_clusters.bed

/home/path/to/software/bedtools/bin/subtractBed -a Intergenic_hg19_all.bed -b intergenic_overlap.bed SV_map.bed > new_intergenics.bed

/home/path/to/software/bedtools/bin/subtractBed -a Clusters_GRCh37.bed -b clusters_overlap.bed SV_map.bed -wao > old_new_clusters.bed

/home/path/to/software/bedtools/bin/subtractBed -a Intergenic_hg19_all.bed -b intergenic_overlap.bed SV_map.bed -wao > old_new_intergenics.bed

#For clusters
CHRS=`seq 1 22`
SETS=`seq 1 5`
for CHR in $CHRS; do
echo "starting chr${CHR}"
 rm chr${CHR}.vcf.gz.tbi
 gunzip chr${CHR}.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chr${CHR}.vcf -b clusters_overlap.bed SV_map.bed > new_chr${CHR}.vcf
bgzip new_chr${CHR}.vcf
tabix -p vcf new_chr${CHR}.vcf.gz


for SET in $SETS; do
echo "starting chr${CHR} population ${SET}"
rm chr${CHR}_${SET}.vcf.gz.tbi
gunzip chr${CHR}_${SET}.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chr${CHR}_${SET}.vcf -b clusters_overlap.bed SV_map.bed > new_chr${CHR}_${SET}.vcf
bgzip new_chr${CHR}_${SET}.vcf
tabix -p vcf new_chr${CHR}_${SET}.vcf.gz

done
done

#test on X chromosome
SETS=`seq 1 5`
rm chrX.vcf.gz.tbi
gunzip chrX.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chrX.vcf -b clusters_overlap.bed SV_map.bed > new_chrX.vcf
bgzip new_chrX.vcf
tabix -p vcf new_chrX.vcf.gz

for SET in $SETS; do
echo "starting population ${SET}"
rm chrX_${SET}.vcf.gz.tbi
gunzip chrX_${SET}.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chrX_${SET}.vcf -b clusters_overlap.bed SV_map.bed > new_chrX_${SET}.vcf
bgzip new_chrX_${SET}.vcf
tabix -p vcf new_chrX_${SET}.vcf.gz
done


#For intergenic regions
CHRS=`seq 1 22`
SETS=`seq 1 5`
for CHR in $CHRS; do
echo "starting chr${CHR}"
rm chr${CHR}_intergenic.vcf.gz.tbi
gunzip chr${CHR}_intergenic.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chr${CHR}_intergenic.vcf -b intergenic_overlap.bed SV_map.bed > new_chr${CHR}_intergenic.vcf
bgzip new_chr${CHR}_intergenic.vcf
tabix -p vcf new_chr${CHR}_intergenic.vcf.gz

for SET in $SETS; do
echo "starting chr${CHR} population ${SET}"
rm chr${CHR}_${SET}_intergenic.vcf.gz.tbi
gunzip chr${CHR}_${SET}_intergenic.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chr${CHR}_${SET}_intergenic.vcf -b intergenic_overlap.bed SV_map.bed > new_chr${CHR}_${SET}_intergenic.vcf
bgzip new_chr${CHR}_${SET}_intergenic.vcf
tabix -p vcf new_chr${CHR}_${SET}_intergenic.vcf.gz
done
done

#test on X chromosome
rm chrX_intergenic.vcf.gz.tbi
gunzip chrX_intergenic.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chrX_intergenic.vcf -b intergenic_overlap.bed SV_map.bed > new_chrX_intergenic.vcf
bgzip new_chrX_intergenic.vcf
tabix -p vcf new_chrX_intergenic.vcf.gz

SETS=`seq 1 5`
for SET in $SETS; do
echo "starting population ${SET}"
rm chrX_${SET}.vcf_intergenic.gz.tbi
gunzip chrX_${SET}_intergenic.vcf.gz
/home/path/to/software/bedtools/bin/subtractBed -header -a chrX_${SET}_intergenic.vcf -b intergenic_overlap.bed SV_map.bed > new_chrX_${SET}_intergenic.vcf
bgzip new_chrX_${SET}_intergenic.vcf
tabix -p vcf new_chrX_${SET}_intergenic.vcf.gz
done

