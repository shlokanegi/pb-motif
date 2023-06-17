----
# PB-Motif

A method for identifying gene/pseudogene rearrangements with PacBio long reads.

----
## Example usages:

Simulated PMS2/PMS2CL rearrangements:

`python pb-motif.py \ `  
`    -i sample_data/pms2_chimera_04/reads.fq.gz \ `  
`    -m motifs/motifs_pms2-pms2cl.p \ `  
`    -o output_dir/ \ `  
`    --skip-qs`

`--skip-qs` is used to bypass quality-score filtering in this instance because the simulated data uses dummy values.

Clinical CAH sample:

`python pb-motif.py \ `  
`    -i sample_data/cyp21a2_normal_01/reads.p \ `  
`    -m motifs/motifs_cyp21a2-a1p.p \ `  
`    -o output_dir/ \ `  
`    -pc 100`

In this case the input data is a pickle of pre-processed reads where motifs have already been identified. These pickles can be saved using the `-p` option. Additionally, `-pc 100` is used so that PB-Motif only outputs plots corresponding to clusters of reads supported by at least 100 reads.
