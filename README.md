# slurm skmer pipeline for quick analysis

## Order of scripts to get finalized output.

*Place your paired fastq files in numbered folders (e.g. folder_1, folder_2, folder_3) for better parallization.*

| Script Title | Use |
|----|----|
| 1_run-bbmap.slurm | runs bbmap_pipeline.sh in ./other_scripts/ (_bbduk.sh_ and _dedupe.sh_ from the BBtools suite).|
| 2_run-consult.slurm | Archaea and Bactrial decontamination using Consult-II. **WoL: Reference Phylogeny for Microbes (bacteria and archaea) (18 Gb - lightweight and robust)** library recommended.|
| 3_run-kraken2.slurm| Human read removal using standard human library as described in Kraken2 manual. |
| 3_run_compile_hist_data.sh | Compiling Necessary Data for RESPECT _(do not submit to slurm scheduler)_ |
| 4_downsample-reads.slurm | Uses "reformat.sh" from BBtools to downsample reads given an estimated genome size and average read length. | 
| 5_run-skmer-query.slurm | Runs *skmer query* command to generate histograms, stats, and mash files for each sample|
| 6_skmer-distance.slurm  | Runs _skmer distance_ to generate distance matrix|

Notes:
- You may have to change the paths where the scripts or libraries have been placed.
- Install kraken2 via anaconda (or export to PATH if installed via github).
