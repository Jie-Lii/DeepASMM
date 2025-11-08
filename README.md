# DeepASMM: Decoding Cis-Regulatory Mechanisms via Quantification of Motif Autonomy and Contextual Synergy

## DeepASMM
**DeepASMM** (**Deep** Learning Driven **A**utonomous and **S**ynergic **M**otif **M**ining Framework) is a neural network–based approach that quantifies both *autonomous functionality* and *sequence context synergy* of motifs through forward-propagation perturbation analysis. 

For additional details, we kindly invite you to refer to the DeepASMM publication:   [*<ins>DeepASMM: Decoding Cis-Regulatory Mechanisms via Quantification of Motif Autonomy and Contextual Synergy</ins>*](https://github.com/Jie-Lii/DeepASMM).  


### Workflow of DeepASMM
We developed a deep learning–based framework called **DeepASMM** to identify motifs with autonomous effects and sequence context synergy in genomic sequences.

DeepASMM consists of four main steps. **First**, deep learning–based genomic sequence prediction models are constructed to capture regulatory information embedded in sequences. One-hot encoding is used to preprocess sequences, enabling the model to effectively learn predictive rules. **Second**, background sequences are selected from true positive samples, ensuring that the sequences used for motif discovery contain real regulatory signals. **Third**, candidate motif localization is performed by scanning background sequences to identify all occurrences of each motif. **Fourth**, motif functionality assessment is conducted: the motif autonomous functionality score quantifies the intrinsic regulatory effect of a motif by embedding it into empty sequences, while the sequence context synergy score measures how the surrounding sequence context influences the motif’s effect.

<div align=center>
<img height="600" src="imgs/fig1.png">
</div>  

### Experimental Data Introduction
In this study, we used the dataset of our previously developed *maize gene expression prediction* model DeepCBA for the experiments ([Wang et al., 2024](https://www.cell.com/plant-communications/fulltext/S2590-3462(24)00293-1)). This dataset includes chromatin interaction and gene expression data of three tissues (shoot, ear, and tassel) of maize (B73).

The *maize chromatin accessibility prediction* task, also involves the data of three tissues: shoot, ear, and tassel ([Peng et al., 2019](https://www.nature.com/articles/s41467-019-10602-5); [Li et al., 2019](https://www.nature.com/articles/s41467-019-10603-4); [Sun et al., 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02063-7)). For each dataset of chromatin accessibility peaks, we extended to 300bp region based on the central locus as positive samples. Negative samples were randomly selected from the maize B73 reference genome with the same number as positive samples, ensuring no overlap with positive regions. All samples were randomly split into training, validation, and test sets at a ratio of *6:2:2*.

For the *human chromatin accessibility prediction* task, we used the dataset reported in the Basset model ([Kelley et al., 2016](https://genome.cshlp.org/content/26/7/990.full)). This dataset contains *2,071,886* sequences of 600bp covering *164* human cell types. In this dataset, *1,930,000* sequences were randomly selected as the training set, *70,000* as the validation set, and *71,886* as the test set.

## Environment  

### CUDA Environment
If you are running this project using GPU, please configure CUDA and cuDNN according to this version.  
  
|  | Version |
|-----:|---------------|
|    CUDA    |    8.0    |
|    cuDNN    |    11.0    |  


### package Environment 
This project is based on Python 3.8.13. The required environment is as follows:  

|                 |    Version  |
|----------------:|-------------|
|    numpy        |    1.19.5   |
|    pandas       |    1.2.4    |
|    tensorflow   |    2.4.0    |
|    tf-keras-vis |    0.8.4    |
|    biopython    |    1.79     |
|    tqdm         |    4.62.3   |  

Some test cases have also been verified to run on *tensorflow 1.15*.

For more required packages, please refer to the [requirements.txt](requirements.txt) file in this project.

## How to Run
1. For detailed instructions, please refer to the [DeepASMM Manual](https://github.com/Jie-Lii/DeepASMM/tree/main/DeepASMM/README.md) in this repository.

2. **Parallel execution:** DeepASMM supports Python-based multi-processing acceleration. Depending on your hardware configuration, up to **10× speedup** can be achieved.

3. **Demo examples** are available at [DeepASMM Demo](https://github.com/Jie-Lii/DeepASMM/tree/main/demo/). 

4. If multiple motif sequence alignments are required, we recommend using our extended tool, [**TOMTOM Parallelization Tool**](https://github.com/Jie-Lii/parallel_tomtom).  
   This tool is built on Python’s multi-processing framework and wraps the MEME Suite TOMTOM module, achieving up to **100× faster alignment** performance.


## Contact
For questions or suggestions, please reach out: [Li_jie@webmail.hzau.edu.cn](mailto:Li_jie@webmail.hzau.edu.cn), [liujianxiao321@webmail.hzau.edu.cn](mailto:liujianxiao321@webmail.hzau.edu.cn).
