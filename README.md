# Translatome_Bcells_glucose
## 1. Biological experiment outline:


#### TREATMENT
- Triplicates (2x15cm dishes = 30 milion cells per replicate), splited 5 days before 1:2
- Starved for 24 hours at 0.5 mM glucose, then 0.5 hours at 20 mM for HIGH glucose while 0.5 mM for LOW glucose

#### FRACTION selection for sequencing
Polysome profiling 13 fractions divided in:
 - Monosomes: 5-6
 - Light polysomes: 7-9
 - Heavy polysomes: 10-13
 - Total RNA for each replicate

## 2. Sequencing Analysis (Shell commands and scripts)

| Raw_reads         | Monosomes | Light Poly | Heavy Poly |
|-------------------|-----------|------------|------------|
| Low_glucose rep1  | 10814918  | 9745454    | 12051978   |
| Low_glucose rep2  | 10578664  | 13997704   | 9749807    |
| Low_glucose rep3  | 13326344  | 11453671   | 7871513    |
| High_glucose rep1 | 12416163  | 10723806   | 8476376    |
| High_glucose rep2 | 10948574  | 11572776   | 8985906    |
| High_glucose rep3 | 10855263  | 8762845    | 13196915   |

 **run_dirty_STAR.sh**
- Mapping: 
STAR 2.5.4b; parameters: --outFilterMultimapNmax 20 --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outStd Log --seedSearchStartLmaxOverLread 0.5 --winAnchorMultimapNmax 36 --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5
 
| Unique mapped     | Monosomes | Light Poly | Heavy Poly | Total RNA |
|-------------------|-----------|------------|------------|-----------|
| Low_glucose rep1  | 61.24%    | 75.83%     | 74.58%     | 79.24%    |
| Low_glucose rep2  | 56.84%    | 76.06%     | 74.69%     | 78.81%    |
| Low_glucose rep3  | 50.12%    | 68.52%     | 73.96%     | 77.56%    |
| High_glucose rep1 | 60.37%    | 73.56%     | 73.78%     | 78.71%    |
| High_glucose rep2 | 58.57%    | 69.51%     | 80.24%     | 79.24%    |
| High_glucose rep3 | 64.97%    | 68.10%     | 77.64%     | 78.40%    |
 
- TPM quantification: 
RSEM (used) and STRINGTIE (just tried)

- Filtering of the files for protein coding genes

## 3.  Statistical analysis (R script)
ask it to Costas
#### Gene ontology Analysis
ask it to Costas

## 4. RNA features analysis (Python script)
ask to costas


