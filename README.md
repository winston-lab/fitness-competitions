
# Flow cytometry fitness competition analysis

## description

An analysis pipeline for diamide flow cytometry fitness competitions with the following major steps:

- Parse raw flow cytometry data in `.fcs` format into human-readable format, and add metadata.
- Perform basic linear thresholding to get cell counts per sample, and do exploratory data visualization.

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

### required files

- `.fcs` format flow cytometry data

## instructions
**0**. Clone this repository.

```bash
git clone https://github.com/winston-lab/fitness-competitions.git
```

**1**. Create and activate the `snakemake_minimal` virtual environment using conda. The virtual environment creation can take a while.

```bash
# navigate into the pipeline directory
cd fitness-competitions

# create the snakemake_minimal environment
conda env create -v -f envs/snakemake_minimal.yaml

# activate the environment
source activate snakemake_minimal

# to deactivate the environment
# source deactivate
```

**2**. For each experiment, create a `.tsv` sample sheet containing the paths to each `.fcs` data file, and their associated metadata. The samplesheet layout should follow the template provided at `samplesheets/samplesheet_template.tsv`.

**3**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` according to your experimental design. Initially, you may not know what cutoffs to set for calling cells YFP or mCherry positive. The pipeline can later be re-run with parameters determined by the data visualization produced by the initial run.

```bash
# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml    # or use your favorite editor
```

**4**. With the `snakemake_minimal` environment activated, do a dry run of the pipeline to see what files will be created.

```bash
./dryrun.sh
```

**5**. The pipeline can be run on the local machine by running `./localrun.sh`. The first time the pipeline is run, conda will create separate virtual environments for some of the jobs to operate in. The pipeline also supports distributed execution on a cluster. On the HMS O2 cluster, which uses the SLURM job scheduler, entering `sbatch slurm_submit.sh` will submit the pipeline as a single job which spawns individual subjobs as necessary. This can be adapted to other job schedulers and clusters by adapting `slurm_submit.sh`, which submits the pipeline to the cluster, `slurm_status.sh`, which handles detection of job status on the cluster, and `cluster.yaml`, which specifies the resource requests for each type of job.


