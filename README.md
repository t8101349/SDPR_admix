# SDPR_admix
A statistical method to calculate PRS in admixed population. SDPR_admix is still in the early stage of testing. If you have encounbtered a problem, please contact zhou1633@purdue.edu.

## Installation

To install SDPR, you need to first download the repo:

```
git clone https://github.com/eldronzhou/SDPR_admix.git
```

You can then compile SDPR by running `make`. If there is run time error that the shared library "libgsl.so.0" not found, you can fix it by typing

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_SDPR_admix_dir/gsl/lib
```

If there is run time error that the shared library "libgsl.so.25" not found, and you install libgsl.so.28, you can fix it by typing

```
cd /opt/conda/envs/test/lib(your_lib_path)
ln -s libgsl.so.28 libgsl.so.25
cd work_dict(your_project_dictionary)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_SDPR_admix_dir/gsl/lib
```

## Quick start

SDPR_admix can be run from the command line. To see the full list of options, please type

```bash
./SDPR_admix -h
```
Below are the required options.

- vcf (required): Path to the phased genotype file in the gzipped vcf format. Not that currently missing genotypes are not supported. If you use Eagle2 and Rfmix2, there shouldn't be missing genotype in the default setting.
- msp (required): Path to path to the directory containing RFMix2 solved local ancestry files.
- pheno (required): path to the phenotype file. The phenotype will be read from the 3rd column of the specified space- or tab-delimited file. There is no header and NA value can be included.
- covar (required): path to the covariate file. Covariates will be reading from the first column. There is no header for the covariate file.
- out (required): Path to the output file containing estimated effect sizes.
- rho (required): Trans-ethnic genetic correlation output by Admix-kit between 0 and 1. Default is 0.8. 

## Running SDPR_admix

```bash
./SDPR_admix -vcf test/chr22_train.vcf.gz -msp test/chr22_train.msp.tsv -pheno test/train.pheno -covar test/covar.txt -out test/res.txt
```
The output file has the following format:

```
1       833068  rs12562034      G       A       0       0
1       843942  rs4040617       A       G       0       0
...
```
where the columns are chromsome, position, variant ID, effect allele, non-effect allele, effect sizes corrsponding to population 0 in RFMix2 file and effect sizes corrsponding to population 1 in the RFMix2 file.

## Derive the PRS

The following command can be used to calculate ancestry-aware PRS for the admixed population.

```bash
./score -vcf test/chr22_train.vcf.gz -msp test/chr22_train.msp.tsv -score test/res.txt -out test/test.profile
```
If you have bgzipped vcf and RFMix2 local ancestry files for chr1-22 with the name `test/chr[1-22]_train.vcf.gz`, then you can use the following command for iterative calculation over all chromsomes:

```bash
./score -vcf test/chr#_train.vcf.gz -msp test/chr#_train.msp.tsv -score test/res.txt -out test/test.profile
```
## Below is a step-by-step guide to the complete process, tailored to your current situation. It includes the required commands and recommended tools to help you navigate from PLINK → VCF → phasing → local ancestry → SDPR_admix.

# ✅ Starting point：The data you have
mydata.bed

mydata.bim

mydata.fam

GWAS statistics .txt (or .ma)

Phenotype (to be analyzed) + covariates (if any)

# 🧱 Step 1: Data Quality Control (QC)

```bash
plink --bfile mydata \
      --maf 0.01 \
      --geno 0.05 \
      --hwe 1e-6 \
      --make-bed \
      --out mydata.qc
```
This will produce the cleaned-up data:

mydata.qc.bed/bim/fam

# 🔄 Step 2: Convert to VCF format

```bash
plink --bfile mydata.qc \
      --recode vcf bgz \
      --out mydata.qc
```

This will produce：

mydata.qc.vcf.gz

# 🧬 Step 3: Phasing
推薦使用 SHAPEIT4

```bash
shapeit4 --input mydata.qc.vcf.gz \
         --map genetic_map_chrXX.txt \
         --region chrXX \
         --output mydata.phased.chrXX.vcf.gz \
         --thread 8
```
If you have multiple chromosomes, you can run it multiple times or write a bash loop.

🧬 步驟四：Local Ancestry 推斷（使用 RFMix v2）
🔹 1. 準備參考族群（如 1000 Genomes）
下載 VCF

分好族群列表（如 EUR.list, AFR.list, AMR.list）

相同 SNP/位置編號與你的目標資料

🔹 2. 轉換格式給 RFMix
使用工具如 vcf2rfmix.py，輸出：

.alleles（haplotype 資料）

.classes（族群標籤）

```bash
python vcf2rfmix.py --vcf mydata.phased.chrXX.vcf.gz \
                    --ref ref_panel.chrXX.vcf.gz \
                    --ref-labels EUR/AFR.list \
                    --out rfmix_input.chrXX
```
🔹 3. 執行 RFMix
```bash
rfmix -f rfmix_input.alleles \
      -r ref_input.alleles \
      -m rfmix_input.classes \
      -g genetic_map_chrXX.txt \
      -o rfmix_output.chrXX
```
RFMix 輸出：

*.msp.tsv（local ancestry information）

🧾 步驟五：整理 SDPR_admix 輸入資料
檔案名稱	說明
chr22.vcf.gz	相相後的 VCF（phased）
chr22.msp.noheader.txt	RFMix 輸出檔（可去掉 header）
train.pheno.txt	表型檔，欄位：FID IID phenotype
covar.tab.txt	共變數檔，欄位：FID IID cov1 cov2 ...
summary.ma	summary statistics for SDPR_admix

▶️ 步驟六：執行 SDPR_admix
bash
複製
編輯
./SDPR_admix \
  -vcf chr22.vcf.gz \
  -msp chr22.msp.noheader.txt \
  -pheno train.pheno.txt \
  -covar covar.tab.txt \
  -rho 0.9 \
  -out result_chr22.txt
執行後會產出：

每個樣本的 PRS 分數

統計評估報告（若指定）
