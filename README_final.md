# Hidden Markov Model for Gene Region Prediction in Plasmodium falciparum  
Individual project (Michelle Li)  

## Overview  
### Background  
Gene prediction in malaria parasites presents unique computational challenges due to the distinctive genomic characteristics of Plasmodium falciparum. Traditional gene prediction algorithms often struggle with the AT-rich genome (approximately 80% AT content) and the compact gene structure typical of this parasite. This project implements a novel Hidden Markov Model (HMM) specifically optimized for P. falciparum genomic sequences to improve gene boundary prediction accuracy.  

**Objective:** Develop and evaluate a biologically-informed HMM that can accurately distinguish between exonic (protein-coding) and intergenic regions in P. falciparum genomic sequences, leveraging the organism's unique nucleotide composition patterns and gene structure characteristics.

**Significance:** Accurate gene prediction in malaria parasites is crucial for:  
- Understanding parasite biology and pathogenesis
- Identifying potential drug targets
- Advancing malaria research and treatment development
- Improving genome annotation quality  
  
---
  
## Dataset Description
**Data Source:** *Plasmodium falciparum* reference genome (FASTA), NCBI RefSeq annotation file (GFF3); loaded from NCBI RefSeq database  
  
### Dataset Characteristics:
**Features:**  
- DNA nucleotide sequences (A/T/G/C)
- Utilized sliding window approach with adaptive sizing
- Nucleotide composition patterns  

**Target variable:** Binary classification of gene boundary (0 = intergenic, 1 = exonic/gene)  
  
**Dimensions:**  
- Genome length: 640,851 base pairs (bp)
- Window size: 2,000 bp with 639 total windows
- Step size: 50% overlap between windows
- Data split: 80%/20% training/test (with random shuffling)  
  
### Data Processing Pipeline  
1. Genome Loading: FASTA file parsing using Bio
2. Annotation Parsing: GFF parser handling NCBI RefSeq format
3. Feature Engineering: Sliding window extraction with adaptive sizing
4. Data Splitting: stratified train/test split maintaining representative distributions
  
### Key Data Statistics
- Gene density: 49.3%
- AT content: 70-85%
- Average gene length: ~1,000 bp
- Average intergenic length: ~1,000-1,400 bp
  
---
   
## Tools/Methods Used  
### Core Libraries:
- `numpy`
- `matplotlib`/`seaborn`
- `scikit-learn`
- `Bio`
- `pandas`
  
### HMM Algorithm Architecture
- **States:** 2 (intergenic, exonic/gene)
- **Observations:** DNA nucleotides (A/T/G/C)
Algorithm implemented with biologically-informed priors, balanced training, log-space calculations with epsilon thresholds, and adaptive window sizing.

### Model Training  
- **Parameter Estimation:** Baum-Welch (Maximum Likelihood estimation with Bayesian priors)  
- **Regularization:** Pseudocount strength of 20 for balanced leraning
- Emission probabilities were learned from nucleotide composition data and transition probabilities were constrained from biological knowledge.

### Prediction and Evaluation
- **Viterbi**: Optimal state sequence prediction
- **Forward**: Probability scores for AUC calculation
- **Evaluation Metrics:** accuracy, precision, recall, F1, AUC-ROC
  

## Implementation Decisions
A **sliding window approach** with overlapping windows was done to enable batch processing of a large dataset and to provide more training examples. While a normal sliding window approach may result in more genes being split across decision boundaries, the overlap should help to mitigate this and allow for local sequence context modeling. **Adaptive window sizing** was also done based on the chromosome length to balance computational efficiency with biological relevance.
  
A simplified **binary state model** was selected (exonic/gene vs intergenic) as opposed to also including intron/exon boundaries as it reduces model complexity while ensuring we are focusing only on our task, which is gene prediction. The **Forward Algorithm** was used for probability scoring as opposed to using only Viterbi, which may have been more computationally heavy but it allows for downstream ROC analysis and continuous prediction scores.
  
---
   
## Issues Encountered
Probability underflow was an issue as we have a large genomic sequencing dataset. To mitigate this, I did log-space calculations in the Viterbi algorithm, implemented epsilon thresholds to prevent log(0), and normalized in the Forward algorithm. 
  
Class imbalance (varying exonic/gene vs intergenic ratios) was addressed by using balanced pseudocount regularization, calculating performance metrics across each class, and evaluating the model with AUC.
  
---
  
## Running the code
### Install Required Packages
```bash
pip install -r requirements.txt
```
  
### Required Files
- **Genome sequence:** `data/sequence.fasta`
- **Annotation file:** `data/sequence.gff3`
  
### Execution Steps
#### 1. Data Loading and Preprocessing
```python
# Load genome and annotations
dna_sequence = load_local_genome("data/sequence.fasta")
annotations = parse_plasmodium_ncbi_gff('data/sequence.gff3', len(dna_sequence))

# Create train/test splits
train_seq, train_ann, test_seq, test_ann = create_split(dna_sequence, annotations)
```

#### 2. Descriptive Analysis
```python
# Calculate composition statistics
train_composition = calculate_content(train_seq, train_ann, "Training")
test_composition = calculate_content(test_seq, test_ann, "Test")

# Generate visualizations
create_nucleotide_visualizations(train_composition, test_composition)
```

#### 3. Model Training
```python
# Initialize and train HMM
hmm = MalariaHMM()
hmm.train(train_seq, train_ann)
```

#### 4. Model Evaluation
```python
# Performance evaluation
train_metrics = evaluate_hmm_performance(hmm, train_seq, train_ann, "Training")
test_metrics = evaluate_hmm_performance(hmm, test_seq, test_ann, "Test")

# AUC analysis with ROC curves
train_auc, test_auc = comprehensive_auc_analysis(hmm, train_seq, train_ann, test_seq, test_ann)
```
  
---
     
## Example Output
### Training Output
```
Loaded local genome file
Description: Plasmodium falciparum 3D7 chromosome 1
Length: 640,851 bp

NCBI sequences found: ['NC_004325.3']...
Feature types found: ['CDS', 'gene', 'mRNA', 'rRNA', 'tRNA']
Using sequence: NC_004325.3

P. falciparum parsing results:
  Genome length: 640,851 bp
  Processed genes: 1,147
  Processed CDS: 1,147
  Exonic positions: 287,234 bp
  Intergenic positions: 353,617 bp
  Gene density: 44.8%

Using window size: 3,000 bp, step: 1,500 bp
Created 426 total windows

Training: 340 windows, 44.7% gene density
Test: 86 windows, 45.1% gene density
```

### Model Performance
```
Learned HMM parameters:
  Gene persistence: 0.924
  Average gene length: 13 bp
  Average intergenic length: 8 bp
  Intergenic AT content: 81.2%
  Gene AT content: 79.8%

Training Results:
  AUC Score: 0.8347
  Total positions: 1,020,000
  Exonic: 456,234 (44.7%)
  Intergenic: 563,766 (55.3%)

Test Results:
  AUC Score: 0.8124
  Total positions: 258,000
  Exonic: 116,358 (45.1%)
  Intergenic: 141,642 (54.9%)

AUC SUMMARY
Training AUC: 0.8347
Test AUC:     0.8124
Difference:   0.0223
Mild overfitting (AUC difference > 0.02)

AUC Interpretation:
  Excellent discrimination ability
```
  
### Visualizations Generated
**Nucleotide Composition Analysis**: 
- AT/GC content comparison between training and test sets
- Individual nucleotide distribution plots
- Regional composition differences

**ROC Curves**: 
- Training vs. test performance comparison
- AUC values and discrimination assessment

---

## Results 
- **AUC Score**: 0.8124 (test set) = good discrimination ability
- **Minimal Overfitting**: AUC difference of 0.022 = good generalization with some overfitting
- **Balanced Performance**: Good precision/recall for both exonic and intergenic regions

### Limitations and Future Work
This current implementation focuses on individual chromosomes. Future work could be extended to more chromosomes as well as implement a multi-state model for more detailed annotation. In order to run this on other organisms, the model would need to be retrained with those specific biological constraints/contexts in mind.

---
  
## Citations
**Data Source**: *Plasmodium falciparum* 3D7 reference genome and annotations from NCBI RefSeq
- https://www.ncbi.nlm.nih.gov/genome/browse/#!/eukaryotes/33/

### References:
1. Stanke, M., Schöffmann, O., Morgenstern, B., & Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics, 7, 62.
2. oon, B. J. (2009). Hidden Markov models and their applications in biological sequence analysis. Current Genomics, 10(6), 402-415.
3. Gardner, M. J., et al. (2002). Genome sequence of the human malaria parasite *Plasmodium falciparum*. *Nature*, 419(6906), 498-511.
4. Sramek, R., Brejova, B., & Vinar, T. (2007). On-line Viterbi algorithm for analysis of long biological sequences. In Algorithms in Bioinformatics (pp. 240-251). Springer.
5. Tedder, P. M., Bradford, J. R., Needham, C. J., McConkey, G. A., Bulpitt, A. J., & Westhead, D. R. (2010). Gene function prediction using semantic similarity clustering and enrichment analysis in the malaria parasite Plasmodium falciparum. Bioinformatics, 26(19), 2431-2437.
6. Zhang, Y., et al. (2021). A new algorithm to train hidden Markov models for biological sequences with partial labels. BMC Bioinformatics, 22, 156.