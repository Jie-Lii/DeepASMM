# Quick Start

## Run with configuration file

```bash
python main.py --txt config.txt
```

## Run with command line arguments

```bash
python main.py
  --input_file_path input.fa
  --model_path model.h5
  --output XXX
  --motif_length 10
```

If `--txt` is provided, all parameters will be loaded from the configuration file and other arguments are not required.


## Demo Examples

You can directly run the following demo datasets:

```bash
python main.py --txt ../demos/basset_demo/command.txt

python main.py --txt ../demos/maize_demo/command.txt
```
* `basset_demo`: Example based on the Basset model. This is a lightweight demo for quickly testing environment configuration, with an approximate runtime of 10 minutes.
* `maize_demo`: Example for maize genomic sequence analysis.


## Arguments

|Parameter|Description|
|-|-|
|`--txt`|Path to parameter configuration file. If provided, all other arguments are optional.|
|`--input_file_path`|Path to input sequence file.|
|`--model_path`|Path to trained model.|
|`--output`|Output directory.|
|`--motif_length`|Length of motifs to mine. Recommended: 8–12 bp.|
|`--num_processor`|Number of parallel processes. Default: 8.|
|`--min_motif_count`|Minimum occurrence count of a motif in background sequences.|
|`--top_r`|Select top r% motifs for quantification (`all` supported).|
|`--min_top_r_num`|Minimum number of motifs retained after filtering.|
|`--category`|Output category used for multi-class models.|
|`--CUDA_device`|CUDA device ID. Set `-1` for CPU.|
|`--chunk_size`|Batch size for inference.|


## Contact

For questions or suggestions, please open an issue in this repository.

