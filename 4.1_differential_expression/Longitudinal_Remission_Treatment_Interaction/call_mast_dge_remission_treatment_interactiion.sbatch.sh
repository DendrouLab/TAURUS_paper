mast_dge_script=/well/cartography/projects/analysis/302_external_oxford_datasets/TAURUS/revisions/src/run.sbatch.mast.dge.treatment.interaction.remission.With.Site.sh
submission_file=/well/cartography/projects/analysis/302_external_oxford_datasets/TAURUS/revisions/src/MAST_DGE_CD_trt_remi_int_with_site_bucket_celltypes_df.tsv
n_channels=`wc -l $submission_file  | tr ' ' '\n' | head -n 1`
# n_channels=4
sbatch --array 2-$n_channels  $mast_dge_script $submission_file
