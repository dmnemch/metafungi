import glob

config = {'threads': 8}
CAT_DB = "../db/CAT_prepare_20210107/2021-01-07_CAT_database/"
CAT_TAX = "../db/CAT_prepare_20210107/2021-01-07_taxonomy/"
KRAKEN2_DB = "../db/kraken2/pluspf"
ITS_DB_BASE = "../db/UNITE/UNITE_public_10.05.2021"
david_wgs_samples=[sample.split("/")[-1].split("_1")[0] for sample in glob.glob(f"../data/david_wgs/raw/*UDB*_1*")] 

rule all:
  input:
    expand("../out/{dataset}/kraken2/{sample}/{sample}.txt", dataset="david_wgs", sample=david_wgs_samples),
    expand("../out/{dataset}/its/{sample}/{sample}_UNITE_public_10.05.2021.sorted.bam", dataset="david_wgs", sample=david_wgs_samples),
    expand("../out/{dataset}/cat/{sample}/{sample}.cat.contig2classification.txt", dataset="david_wgs", sample=david_wgs_samples)

# =================== CUTADAPT ===================
rule cutadapt:
  input:
    r1="../data/{dataset}/raw/{sample}_1.fastq.gz",
    r2="../data/{dataset}/raw/{sample}_2.fastq.gz"
  output:
    r1="../data/{dataset}/trimmed/{sample}_1.cutadapt.fq.gz",
    r2="../data/{dataset}/trimmed/{sample}_2.cutadapt.fq.gz"
  conda: "../envs/cutadapt_env.yml"
  threads: 8
  shell: "cutadapt --quality-cutoff 15,15 --minimum-length 1 -o {output.r1} -p {output.r2} --cores={threads} {input.r1} {input.r2}"

# ==================== CATBAT ====================
rule metaspades:
  input:
    r1="../data/{dataset}/trimmed/{sample}_1.cutadapt.fq.gz",
    r2="../data/{dataset}/trimmed/{sample}_2.cutadapt.fq.gz"
  output: "../out/{dataset}/metaspades/{sample}/contigs.fasta"
  conda: "../envs/metaspades_env.yml"
  threads: config['threads']
  params: out_dir="../out/{dataset}/metaspades/{sample}/"
  shell: "spades.py --meta --threads {threads} -1 {input.r1} -2 {input.r2} -o {params.out_dir}"

rule cat:
  input:
    contigs="../out/{dataset}/metaspades/{sample}/contigs.fasta",
    db=CAT_DB,
    tax=CAT_TAX
  output: "../out/{dataset}/cat/{sample}/{sample}.cat.contig2classification.txt"
  params: prefix="../out/{dataset}/cat/{sample}/{sample}.cat"
  threads: config['threads']
  conda: "../envs/catbat_env.yml"
  shell: "CAT contigs -c {input.contigs} -d {input.db} -t {input.tax} -o {params.prefix} --nproc {threads}"

# ==================== KRAKEN2 ===================
rule kraken2:
  input:
    r1="../data/{dataset}/trimmed/{sample}_1.cutadapt.fq.gz",
    r2="../data/{dataset}/trimmed/{sample}_2.cutadapt.fq.gz",
    db=KRAKEN2_DB
  output:
    out="../out/{dataset}/kraken2/{sample}/{sample}.txt",
    report="../out/{dataset}/kraken2/{sample}/{sample}_report.txt"
  threads: config['threads']
  conda: "../envs/kraken2_env.yml"
  shell: "kraken2 --db {input.db} --paired --threads {threads} --use-names --report {output.report} --output {output.out} {input.r1} {input.r2}"

# ====================== ITS =====================
rule its_bowtie2_build:
  input: "../db/UNITE/{ref}.fasta"
  output:
    "../db/UNITE/{ref}.1.bt2",
    "../db/UNITE/{ref}.2.bt2",
    "../db/UNITE/{ref}.3.bt2",
    "../db/UNITE/{ref}.4.bt2",
    "../db/UNITE/{ref}.rev.1.bt2",
    "../db/UNITE/{ref}.rev.2.bt2"
  conda: "../envs/bt2_env.yml"
  threads: config['threads']
  params: base="../db/UNITE/{ref}"
  shell: "bowtie2-build {input} {params.base} --threads {threads}"

rule its_bowtie2:
  input:
    "../db/UNITE/{ref}.1.bt2",
    "../db/UNITE/{ref}.2.bt2",
    "../db/UNITE/{ref}.3.bt2",
    "../db/UNITE/{ref}.4.bt2",
    "../db/UNITE/{ref}.rev.1.bt2",
    "../db/UNITE/{ref}.rev.2.bt2",
    r1="../data/{dataset}/trimmed/{sample}_1.cutadapt.fq.gz",
    r2="../data/{dataset}/trimmed/{sample}_2.cutadapt.fq.gz"
  output: "../out/{dataset}/its/{sample}/{sample}_{ref}.sorted.bam"
  params: base="../db/UNITE/{ref}"
  conda: "../envs/bt2_samtools_env.yml"
  threads: config['threads']
  shell:
    """
    bowtie2 -q --threads {threads} -x {params.base} -1 {input.r1} -2 {input.r2} | samtools view -bS -q 30 - | samtools sort -o {output} -
    samtools index {output}
    """
