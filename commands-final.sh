conda activate qiime2-amplicon-2023.9

cd /Users/zacharyburcham/Dropbox/UTK/Research/Collab_projects/jj_colombia/qiime2-2023.9

mkdir import feature_tables rep_seqs taxonomy metadata rarefaction tree

#importing paired data
qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path import/manifest-pe.txt \
	--output-path import/paired-end-demux.qza \
	--input-format PairedEndFastqManifestPhred33V2

# remove primers
qiime cutadapt trim-paired \
	--i-demultiplexed-sequences import/paired-end-demux.qza \
	--p-cores 4 \
	--p-front-f CCTACGGGNGGCWGCAG \
	--p-front-r GACTACHVGGGTATCTAATCC \
	--p-match-read-wildcards \
	--p-match-adapter-wildcards \
	--o-trimmed-sequences import/paired-end-demux-trimmed.qza
		
#summarize demux data
qiime demux summarize \
  --i-data import/paired-end-demux-trimmed.qza \
  --o-visualization import/paired-end-demux-trimmed.qzv
  
#dada2 denoising
qiime dada2 denoise-paired \
	--i-demultiplexed-seqs import/paired-end-demux-trimmed.qza \
	--p-trunc-len-f 259 \
	--p-trunc-len-r 235 \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-n-threads 4 \
	--o-table feature_tables/feature_table.qza \
	--o-representative-sequences rep_seqs/rep_seqs.qza \
	--o-denoising-stats rep_seqs/dada2-denoise-stats.qza

qiime feature-table summarize \
  --i-table feature_tables/feature_table.qza \
  --m-sample-metadata-file metadata/metadata.tsv \
  --o-visualization feature_tables/feature_table.qzv
  
#denoising stats
qiime metadata tabulate \
  --m-input-file rep_seqs/dada2-denoise-stats.qza  \
  --o-visualization rep_seqs/dada2-denoise-stats.qzv
  
qiime metadata tabulate \
  --m-input-file rep_seqs/rep_seqs.qza  \
  --o-visualization rep_seqs/rep_seqs.qzv
  
# PREBUILT FULL LENGTH SILVA 138 WITH WEIGHTS classify rep seqs
qiime feature-classifier classify-sklearn \
  --i-reads rep_seqs/rep_seqs.qza \
  --i-classifier taxonomy/silva-138-99-nb-weighted-classifier.qza \
  --o-classification taxonomy/taxonomy-weighted.qza \
  --p-n-jobs 4

qiime metadata tabulate \
  --m-input-file taxonomy/taxonomy-weighted.qza \
  --o-visualization taxonomy/taxonomy-weighted.qzv

### fragment insertion phylogenetic tree creation
qiime fragment-insertion sepp \
  --i-representative-sequences rep_seqs/rep_seqs.qza \
  --i-reference-database tree/sepp-refs-silva-128.qza \
  --o-tree tree/tree.qza \
  --o-placements tree/tree_placements.qza \
  --p-threads 4
  
### Filter feature table of mitochondria and chloroplast and keep only bacteria
qiime taxa filter-table \
 --i-table feature_tables/feature_table.qza \
 --i-taxonomy taxonomy/taxonomy-weighted.qza \
 --p-include bacteria \
 --p-exclude mitochondria,chloroplast \
 --o-filtered-table feature_tables/feature_table-no-chlo-mito.qza

### Filter features based on frequency to remove features less than 10 times present and not in at least 2 samples
qiime feature-table filter-features \
 --i-table feature_tables/feature_table-no-chlo-mito.qza \
 --p-min-frequency 10 \
 --p-min-samples 2 \
 --o-filtered-table feature_tables/feature_table-no-chlo-mito-filt.qza

qiime feature-table summarize \
  --i-table feature_tables/feature_table-no-chlo-mito-filt.qza \
  --o-visualization feature_tables/feature_table-no-chlo-mito-filt.qzv \
  --m-sample-metadata-file metadata/metadata.tsv
  
### Filter out replicate samples with the smaller number of seqs
qiime feature-table filter-samples \
  --i-table feature_tables/feature_table-no-chlo-mito-filt.qza \
  --m-metadata-file samples-to-keep.tsv \
  --o-filtered-table feature_tables/final-table.qza

qiime feature-table summarize \
  --i-table feature_tables/final-table.qza \
  --o-visualization feature_tables/final-table.qzv \
  --m-sample-metadata-file metadata/metadata.tsv

### Create taxa barplot
qiime taxa barplot \
  --i-table feature_tables/final-table.qza \
  --i-taxonomy taxonomy/taxonomy-weighted.qza \
  --m-metadata-file metadata/metadata.tsv \
  --o-visualization taxonomy/taxa-bar-plots.qzv
 
### Check diversity metrics with different rarefying depths - so issues so to keep all samples 28,000 should be good! 
qiime diversity alpha-rarefaction \
	--i-table feature_tables/final-table.qza \
	--i-phylogeny tree/tree.qza \
	--p-max-depth 100000 \
	--p-min-depth 1 \
	--p-steps 21 \
	--o-visualization rarefaction/alpha-rarefaction-100000.qzv
	
#### Core metrics at 25000 rarefaction
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny tree/tree.qza \
  --i-table feature_tables/final-table.qza \
  --p-sampling-depth 25000 \
  --m-metadata-file metadata/metadata.tsv \
  --output-dir core-metrics-results-25000 

mkdir core-metrics-results-25000/alpha_diversity
mkdir core-metrics-results-25000/alpha_diversity/group_signif
mkdir core-metrics-results-25000/alpha_diversity/correlation
mkdir core-metrics-results-25000/beta_diversity
mkdir core-metrics-results-25000/beta_diversity/group_signif
mkdir core-metrics-results-25000/beta_diversity/mantel

mv core-metrics-results-25000/*_vector.qza core-metrics-results-25000/alpha_diversity
mv core-metrics-results-25000/*_emperor.qzv core-metrics-results-25000/beta_diversity
mv core-metrics-results-25000/*_pcoa*.qza core-metrics-results-25000/beta_diversity
mv core-metrics-results-25000/*_matrix.qza core-metrics-results-25000/beta_diversity

## Alpha group significance
# Faith's pd
qiime diversity alpha-group-significance \
  --i-alpha-diversity  core-metrics-results-25000/alpha_diversity/faith_pd_vector.qza \
  --m-metadata-file metadata/metadata.tsv \
  --o-visualization  core-metrics-results-25000/alpha_diversity/group_signif/faith-pd-group-significance.qzv

# Peilou's evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity  core-metrics-results-25000/alpha_diversity/evenness_vector.qza \
  --m-metadata-file metadata/metadata.tsv \
  --o-visualization  core-metrics-results-25000/alpha_diversity/group_signif/evenness-significance.qzv

# Shannon
qiime diversity alpha-group-significance \
  --i-alpha-diversity  core-metrics-results-25000/alpha_diversity/shannon_vector.qza \
  --m-metadata-file metadata/metadata.tsv \
  --o-visualization  core-metrics-results-25000/alpha_diversity/group_signif/shannon-significance.qzv

# Observed features (richness)
qiime diversity alpha-group-significance \
  --i-alpha-diversity  core-metrics-results-25000/alpha_diversity/observed_features_vector.qza \
  --m-metadata-file metadata/metadata.tsv \
  --o-visualization  core-metrics-results-25000/alpha_diversity/group_signif/observed_otus_significance.qzv

## Beta Diversity Group Significance
## Unweighted
# Infection status
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-25000/beta_diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-column infection_status \
  --o-visualization core-metrics-results-25000/beta_diversity/group_signif/infection_status-unweighted-unifrac-significance.qzv 

# Symptomatic
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-25000/beta_diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-column symptomatic \
  --o-visualization core-metrics-results-25000/beta_diversity/group_signif/symptomatic-unweighted-unifrac-significance.qzv 
  
# Group
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-25000/beta_diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-column group \
  --o-visualization core-metrics-results-25000/beta_diversity/group_signif/group-unweighted-unifrac-significance.qzv \
  --p-pairwise

# adonis testing
mkdir core-metrics-results-25000/beta_diversity/adonis

# all metadata groups
qiime diversity adonis \
	--i-distance-matrix core-metrics-results-25000/beta_diversity/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status+symptomatic+group" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization core-metrics-results-25000/beta_diversity/adonis/unweighted-all-adonis.qzv
	
# infection status and symptomatic 
qiime diversity adonis \
	--i-distance-matrix core-metrics-results-25000/beta_diversity/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status*symptomatic" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization core-metrics-results-25000/beta_diversity/adonis/unweighted-infection_status-symptomatic-adonis.qzv
	
## Weighted 
# Infection status
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-25000/beta_diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-column infection_status \
  --o-visualization core-metrics-results-25000/beta_diversity/group_signif/infection_status-weighted-unifrac-significance.qzv 

# Symptomatic
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-25000/beta_diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-column symptomatic \
  --o-visualization core-metrics-results-25000/beta_diversity/group_signif/symptomatic-weighted-unifrac-significance.qzv 
  
# Group
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-25000/beta_diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-column group \
  --o-visualization core-metrics-results-25000/beta_diversity/group_signif/group-weighted-unifrac-significance.qzv \
  --p-pairwise
  
# adonis testing
# all metadata groups
qiime diversity adonis \
	--i-distance-matrix core-metrics-results-25000/beta_diversity/weighted_unifrac_distance_matrix.qza \
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status+symptomatic+group" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization core-metrics-results-25000/beta_diversity/adonis/weighted-all-adonis.qzv
	
# infection status and symptomatic 
qiime diversity adonis \
	--i-distance-matrix core-metrics-results-25000/beta_diversity/weighted_unifrac_distance_matrix.qza \
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status*symptomatic" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization core-metrics-results-25000/beta_diversity/adonis/weighted-infection_status-symptomatic-adonis.qzv
	
### Collapse tables by taxonomy
# Genus
qiime taxa collapse \
 --i-table feature_tables/final-table.qza \
 --i-taxonomy taxonomy/taxonomy-weighted.qza \
 --p-level 6 \
 --o-collapsed-table feature_tables/genus-final-table.qza

qiime feature-table summarize \
  --i-table feature_tables/genus-final-table.qza \
  --o-visualization feature_tables/genus-final-table.qzv \
  --m-sample-metadata-file metadata/metadata.tsv
  
# Family
qiime taxa collapse \
 --i-table feature_tables/final-table.qza \
 --i-taxonomy taxonomy/taxonomy-weighted.qza \
 --p-level 5 \
 --o-collapsed-table feature_tables/family-final-table.qza
 
qiime feature-table summarize \
  --i-table feature_tables/family-final-table.qza \
  --o-visualization feature_tables/family-final-table.qzv \
  --m-sample-metadata-file metadata/metadata.tsv
  
# Order
qiime taxa collapse \
 --i-table feature_tables/final-table.qza \
 --i-taxonomy taxonomy/taxonomy-weighted.qza \
 --p-level 4 \
 --o-collapsed-table feature_tables/order-final-table.qza
 
qiime feature-table summarize \
  --i-table feature_tables/order-final-table.qza \
  --o-visualization feature_tables/order-final-table.qzv \
  --m-sample-metadata-file metadata/metadata.tsv
  
# Class
qiime taxa collapse \
 --i-table feature_tables/final-table.qza \
 --i-taxonomy taxonomy/taxonomy-weighted.qza \
 --p-level 3 \
 --o-collapsed-table feature_tables/class-final-table.qza
 
qiime feature-table summarize \
  --i-table feature_tables/class-final-table.qza \
  --o-visualization feature_tables/class-final-table.qzv \
  --m-sample-metadata-file metadata/metadata.tsv
  
# Phylum
qiime taxa collapse \
 --i-table feature_tables/final-table.qza \
 --i-taxonomy taxonomy/taxonomy-weighted.qza \
 --p-level 2 \
 --o-collapsed-table feature_tables/phylum-final-table.qza
 
qiime feature-table summarize \
  --i-table feature_tables/phylum-final-table.qza \
  --o-visualization feature_tables/phylum-final-table.qzv \
  --m-sample-metadata-file metadata/metadata.tsv
  
### ANCOM-BC
mkdir ancombc ancombc/single ancombc/single/group ancombc/single/infection_status ancombc/single/symptomatic ancombc/multi/

## single formula
## group
# ASV
qiime composition ancombc \
    --i-table feature_tables/final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/single/group/asv-group-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/group/asv-group-diff.qza \
	--o-visualization ancombc/single/group/asv-group-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/group/asv-group-diff.qza \
	--o-visualization ancombc/single/group/asv-group-da-barplot.qzv

# genus
qiime composition ancombc \
    --i-table feature_tables/genus-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/single/group/genus-group-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/group/genus-group-diff.qza \
	--o-visualization ancombc/single/group/genus-group-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/group/genus-group-diff.qza \
	--o-visualization ancombc/single/group/genus-group-da-barplot.qzv
	 
# family
qiime composition ancombc \
    --i-table feature_tables/family-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/single/group/family-group-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/group/family-group-diff.qza \
	--o-visualization ancombc/single/group/family-group-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/group/family-group-diff.qza \
	--o-visualization ancombc/single/group/family-group-da-barplot.qzv
	
# order
qiime composition ancombc \
    --i-table feature_tables/order-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/single/group/order-group-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/group/order-group-diff.qza \
	--o-visualization ancombc/single/group/order-group-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/group/order-group-diff.qza \
	--o-visualization ancombc/single/group/order-group-da-barplot.qzv
	
# class
qiime composition ancombc \
    --i-table feature_tables/class-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/single/group/class-group-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/group/class-group-diff.qza \
	--o-visualization ancombc/single/group/class-group-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/group/class-group-diff.qza \
	--o-visualization ancombc/single/group/class-group-da-barplot.qzv
	
# phylum
qiime composition ancombc \
    --i-table feature_tables/phylum-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/single/group/phylum-group-diff.qza
  
qiime composition tabulate \
	--i-data ancombc/single/group/phylum-group-diff.qza \
	--o-visualization ancombc/single/group/phylum-group-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/group/phylum-group-diff.qza \
	--o-visualization ancombc/single/group/phylum-group-da-barplot.qzv
	
## infection_status
# ASV
qiime composition ancombc \
    --i-table feature_tables/final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/single/infection_status/asv-infection_status-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/infection_status/asv-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/asv-infection_status-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/infection_status/asv-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/asv-infection_status-da-barplot.qzv

# genus
qiime composition ancombc \
    --i-table feature_tables/genus-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/single/infection_status/genus-infection_status-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/infection_status/genus-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/genus-infection_status-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/infection_status/genus-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/genus-infection_status-da-barplot.qzv
	 
# family
qiime composition ancombc \
    --i-table feature_tables/family-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/single/infection_status/family-infection_status-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/infection_status/family-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/family-infection_status-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/infection_status/family-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/family-infection_status-da-barplot.qzv
	
# order
qiime composition ancombc \
    --i-table feature_tables/order-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/single/infection_status/order-infection_status-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/infection_status/order-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/order-infection_status-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/infection_status/order-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/order-infection_status-da-barplot.qzv
	
# class
qiime composition ancombc \
    --i-table feature_tables/class-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/single/infection_status/class-infection_status-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/infection_status/class-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/class-infection_status-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/infection_status/class-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/class-infection_status-da-barplot.qzv
	
# phylum
qiime composition ancombc \
    --i-table feature_tables/phylum-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/single/infection_status/phylum-infection_status-diff.qza
  
qiime composition tabulate \
	--i-data ancombc/single/infection_status/phylum-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/phylum-infection_status-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/infection_status/phylum-infection_status-diff.qza \
	--o-visualization ancombc/single/infection_status/phylum-infection_status-da-barplot.qzv
  
## symptomatic
# ASV
qiime composition ancombc \
    --i-table feature_tables/final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
     --p-reference-levels "symptomatic::no" \
     --p-conserve \
    --o-differentials ancombc/single/symptomatic/asv-symptomatic-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/symptomatic/asv-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/asv-symptomatic-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/symptomatic/asv-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/asv-symptomatic-da-barplot.qzv

# genus
qiime composition ancombc \
    --i-table feature_tables/genus-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
     --p-reference-levels "symptomatic::no" \
     --p-conserve \
    --o-differentials ancombc/single/symptomatic/genus-symptomatic-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/symptomatic/genus-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/genus-symptomatic-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/symptomatic/genus-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/genus-symptomatic-da-barplot.qzv
	 
# family
qiime composition ancombc \
    --i-table feature_tables/family-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
     --p-reference-levels "symptomatic::no" \
     --p-conserve \
    --o-differentials ancombc/single/symptomatic/family-symptomatic-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/symptomatic/family-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/family-symptomatic-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/symptomatic/family-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/family-symptomatic-da-barplot.qzv
	
# order
qiime composition ancombc \
    --i-table feature_tables/order-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
     --p-reference-levels "symptomatic::no" \
     --p-conserve \
    --o-differentials ancombc/single/symptomatic/order-symptomatic-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/symptomatic/order-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/order-symptomatic-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/symptomatic/order-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/order-symptomatic-da-barplot.qzv
	
# class
qiime composition ancombc \
    --i-table feature_tables/class-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
     --p-reference-levels "symptomatic::no" \
     --p-conserve \
    --o-differentials ancombc/single/symptomatic/class-symptomatic-diff.qza

qiime composition tabulate \
	--i-data ancombc/single/symptomatic/class-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/class-symptomatic-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/symptomatic/class-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/class-symptomatic-da-barplot.qzv
	
# phylum
qiime composition ancombc \
    --i-table feature_tables/phylum-final-table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
    --p-reference-levels "symptomatic::no" \
    --p-conserve \
    --o-differentials ancombc/single/symptomatic/phylum-symptomatic-diff.qza
  
qiime composition tabulate \
	--i-data ancombc/single/symptomatic/phylum-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/phylum-symptomatic-diff-tabulate.qzv
	
qiime composition da-barplot \
	--i-data ancombc/single/symptomatic/phylum-symptomatic-diff.qza \
	--o-visualization ancombc/single/symptomatic/phylum-symptomatic-da-barplot.qzv
   
### create relative abundance tables
# phylum
qiime feature-table relative-frequency \
	--i-table feature_tables/phylum-final-table.qza \
	--o-relative-frequency-table feature_tables/phylum-final-table-relfreq.qza
	
qiime feature-table transpose \
	--i-table feature_tables/phylum-final-table-relfreq.qza \
	--o-transposed-feature-table feature_tables/phylum-final-table-relfreq-trans.qza
	
qiime tools export \
  --input-path feature_tables/phylum-final-table-relfreq-trans.qza \
  --output-path feature_tables/phylum-final-table-relfreq-trans
  
mv feature_tables/phylum-final-table-relfreq-trans/feature-table.biom feature_tables/phylum-final-table-relfreq-trans.biom
rmdir feature_tables/phylum-final-table-relfreq-trans/
biom convert -i feature_tables/phylum-final-table-relfreq-trans.biom -o feature_tables/phylum-final-table-relfreq-trans.tsv --to-tsv

# class
qiime feature-table relative-frequency \
	--i-table feature_tables/class-final-table.qza \
	--o-relative-frequency-table feature_tables/class-final-table-relfreq.qza
	
qiime feature-table transpose \
	--i-table feature_tables/class-final-table-relfreq.qza \
	--o-transposed-feature-table feature_tables/class-final-table-relfreq-trans.qza
	
qiime tools export \
  --input-path feature_tables/class-final-table-relfreq-trans.qza \
  --output-path feature_tables/class-final-table-relfreq-trans
  
mv feature_tables/class-final-table-relfreq-trans/feature-table.biom feature_tables/class-final-table-relfreq-trans.biom
rmdir feature_tables/class-final-table-relfreq-trans/
biom convert -i feature_tables/class-final-table-relfreq-trans.biom -o feature_tables/class-final-table-relfreq-trans.tsv --to-tsv

# order
qiime feature-table relative-frequency \
	--i-table feature_tables/order-final-table.qza \
	--o-relative-frequency-table feature_tables/order-final-table-relfreq.qza
	
qiime feature-table transpose \
	--i-table feature_tables/order-final-table-relfreq.qza \
	--o-transposed-feature-table feature_tables/order-final-table-relfreq-trans.qza
	
qiime tools export \
  --input-path feature_tables/order-final-table-relfreq-trans.qza \
  --output-path feature_tables/order-final-table-relfreq-trans
  
mv feature_tables/order-final-table-relfreq-trans/feature-table.biom feature_tables/order-final-table-relfreq-trans.biom
rmdir feature_tables/order-final-table-relfreq-trans/
biom convert -i feature_tables/order-final-table-relfreq-trans.biom -o feature_tables/order-final-table-relfreq-trans.tsv --to-tsv


###################################################### metabolite diversity
# import metabolite data
biom convert -i metabolites/20171207_table.tsv -o metabolites/20171207_table.biom --table-type="OTU table" --to-hdf5
biom convert -i metabolites/180328_table.tsv -o metabolites/180328_table.biom --table-type="OTU table" --to-hdf5

qiime tools import \
        --input-path metabolites/20171207_table.biom \
        --output-path metabolites/20171207_table.qza \
        --type 'FeatureTable[Frequency]'

qiime tools import \
        --input-path metabolites/180328_table.biom \
        --output-path metabolites/180328_table.qza \
        --type 'FeatureTable[Frequency]'        

# merge the two runs, rename to match microbial samples names
qiime feature-table merge \
	--i-tables metabolites/20171207_table.qza metabolites/180328_table.qza \
	--o-merged-table metabolites/full_metabolite_table.qza

# remove the CCN11A replicate and keep CCN11B due to higher total abundance
qiime feature-table filter-samples \
  --i-table metabolites/full_metabolite_table.qza \
  --m-metadata-file metadata/metab_metadata.tsv \
  --p-where '[SampleID]="CCN11A"' \
  --p-exclude-ids \
  --o-filtered-table metabolites/full_metabolite_norep_table.qza
  
# rename to match microbial samples names
qiime feature-table rename-ids \
	--i-table metabolites/full_metabolite_norep_table.qza \
	--m-metadata-file metabolites/rename.txt \
	--m-metadata-column new_name \
	--o-renamed-table metabolites/full_metabolite_rename_table.qza

# remove metabolites that were not present in both runs
qiime feature-table filter-features \
  --i-table metabolites/full_metabolite_rename_table.qza \
  --m-metadata-file metabolites/shared-features.tsv \
  --o-filtered-table metabolites/full_metabolite_shared_table.qza

# summarize table
qiime feature-table summarize \
	--i-table metabolites/full_metabolite_shared_table.qza \
	--o-visualization metabolites/full_metabolite_shared_table.qzv

###### export metabolite table, run scaling
qiime tools export \
  --input-path metabolites/full_metabolite_shared_table.qza \
  --output-path metabolites/full_metabolite_shared_table

mv metabolites/full_metabolite_shared_table/feature-table.biom metabolites/full_metabolite_shared_table.biom
rmdir metabolites/full_metabolite_shared_table/
biom convert -i metabolites/full_metabolite_shared_table.biom -o metabolites/full_metabolite_shared_table.tsv --to-tsv

# reimport new table
biom convert -i metabolites/full_metabolite_shared_scaled_table.tsv -o metabolites/full_metabolite_shared_scaled_table.biom --table-type="OTU table" --to-hdf5

qiime tools import \
        --input-path metabolites/full_metabolite_shared_scaled_table.biom \
        --output-path metabolites/full_metabolite_shared_scaled_table.qza \
        --type 'FeatureTable[Frequency]'

# summarize table
qiime feature-table summarize \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--o-visualization metabolites/full_metabolite_shared_scaled_table.qzv

# diversity metrics
mkdir metabolites/diversity metabolites/diversity/alpha metabolites/diversity/beta

# alpha div calculation
qiime diversity alpha \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'observed_features' \
	--o-alpha-diversity metabolites/diversity/alpha/observed_features_vector.qza
	
qiime diversity alpha \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'pielou_e' \
	--o-alpha-diversity metabolites/diversity/alpha/evenness_vector.qza
	
qiime diversity alpha \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'shannon' \
	--o-alpha-diversity metabolites/diversity/alpha/shannon_vector.qza

# alpha group significance
# Peilou's evenness
qiime diversity alpha-group-significance \
	--i-alpha-diversity  metabolites/diversity/alpha/evenness_vector.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization  metabolites/diversity/alpha/evenness-significance.qzv

# Shannon
qiime diversity alpha-group-significance \
	--i-alpha-diversity  metabolites/diversity/alpha/shannon_vector.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization  metabolites/diversity/alpha/shannon-significance.qzv

# Observed features (richness)
qiime diversity alpha-group-significance \
	--i-alpha-diversity  metabolites/diversity/alpha/observed_features_vector.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization  metabolites/diversity/alpha/observed_otus_significance.qzv

# beta div distances
qiime diversity beta \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'braycurtis' \
	--o-distance-matrix metabolites/diversity/beta/bc-dist.qza
	
qiime diversity beta \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'jaccard' \
	--o-distance-matrix metabolites/diversity/beta/jc-dist.qza
	
qiime diversity beta \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'canberra' \
	--o-distance-matrix metabolites/diversity/beta/canberra-dist.qza
	
qiime diversity beta \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--p-metric 'aitchison' \
	--o-distance-matrix metabolites/diversity/beta/aitchison-dist.qza

# pcoa results
qiime diversity pcoa \
	--i-distance-matrix metabolites/diversity/beta/bc-dist.qza \
	--o-pcoa metabolites/diversity/beta/bc-pcoa.qza
	
qiime diversity pcoa \
	--i-distance-matrix metabolites/diversity/beta/jc-dist.qza \
	--o-pcoa metabolites/diversity/beta/jc-pcoa.qza
	
qiime diversity pcoa \
	--i-distance-matrix metabolites/diversity/beta/canberra-dist.qza \
	--o-pcoa metabolites/diversity/beta/canberra-pcoa.qza
	
qiime diversity pcoa \
	--i-distance-matrix metabolites/diversity/beta/aitchison-dist.qza \
	--o-pcoa metabolites/diversity/beta/aitchison-pcoa.qza
	
# apply tnse to canberra
qiime diversity tsne \
	--i-distance-matrix metabolites/diversity/beta/canberra-dist.qza \
	--p-random-state 42 \
	--o-tsne metabolites/diversity/beta/canberra-tnse-pcoa.qza
	
# apply umap to canberra
qiime diversity umap \
	--i-distance-matrix metabolites/diversity/beta/canberra-dist.qza \
	--p-random-state 42 \
	--o-umap metabolites/diversity/beta/canberra-umap-pcoa.qza
	
# apply tnse to aitchison
qiime diversity tsne \
	--i-distance-matrix metabolites/diversity/beta/aitchison-dist.qza \
	--p-random-state 42 \
	--o-tsne metabolites/diversity/beta/aitchison-tnse-pcoa.qza
	
# apply umap to aitchison
qiime diversity umap \
	--i-distance-matrix metabolites/diversity/beta/aitchison-dist.qza \
	--p-random-state 42 \
	--o-umap metabolites/diversity/beta/aitchison-umap-pcoa.qza

# emperor plots
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/bc-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/bc-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/jc-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/jc-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/canberra-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/canberra-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/canberra-tnse-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/canberra-tnse-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/canberra-umap-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/canberra-umap-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/aitchison-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/aitchison-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/aitchison-tnse-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/aitchison-tnse-emperor.qzv
	
qiime emperor plot \
	--i-pcoa metabolites/diversity/beta/aitchison-umap-pcoa.qza \
	--m-metadata-file metadata/metadata.tsv \
	--o-visualization metabolites/diversity/beta/aitchison-umap-emperor.qzv

# adonis testing
mkdir metabolites/diversity/beta/adonis/

# all metadata groups
# bray-curtis
qiime diversity adonis \
	--i-distance-matrix metabolites/diversity/beta/bc-dist.qza\
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status+symptomatic+group" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization metabolites/diversity/beta/adonis/bc-all-adonis.qzv

# jaccard
qiime diversity adonis \
	--i-distance-matrix metabolites/diversity/beta/jc-dist.qza\
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status+symptomatic+group" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization metabolites/diversity/beta/adonis/jc-all-adonis.qzv

# canberra
qiime diversity adonis \
	--i-distance-matrix metabolites/diversity/beta/canberra-dist.qza\
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status+symptomatic+group" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization metabolites/diversity/beta/adonis/canberra-all-adonis.qzv

# aitchison
qiime diversity adonis \
	--i-distance-matrix metabolites/diversity/beta/aitchison-dist.qza\
	--m-metadata-file metadata/metadata.tsv \
	--p-formula "infection_status+symptomatic+group" \
	--p-permutations 999 \
	--p-n-jobs 4 \
	--o-visualization metabolites/diversity/beta/adonis/aitchison-all-adonis.qzv
	
# Metabolites ANCOM-BC
### ANCOM-BC
mkdir ancombc/metabolites
  
## group
qiime composition ancombc \
    --i-table metabolites/full_metabolite_shared_scaled_table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula group \
    --p-reference-levels "group::CUA" \
    --p-conserve \
    --o-differentials ancombc/metabolites/group-diff.qza

qiime composition tabulate \
	--i-data ancombc/metabolites/group-diff.qza \
	--o-visualization ancombc/metabolites/group-diff.qzv
	
qiime composition da-barplot \
	--i-data ancombc/metabolites/group-diff.qza \
	--o-visualization ancombc/metabolites/group-da-barplot.qzv
	
## infection_status
qiime composition ancombc \
    --i-table metabolites/full_metabolite_shared_scaled_table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula infection_status \
    --p-reference-levels "infection_status::uninfected" \
    --p-conserve \
    --o-differentials ancombc/metabolites/infection_status-diff.qza

qiime composition tabulate \
	--i-data ancombc/metabolites/infection_status-diff.qza \
	--o-visualization ancombc/metabolites/infection_status-diff.qzv
	
qiime composition da-barplot \
	--i-data ancombc/metabolites/infection_status-diff.qza \
	--o-visualization ancombc/metabolites/infection_status-da-barplot.qzv

## symptomatic
qiime composition ancombc \
    --i-table metabolites/full_metabolite_shared_scaled_table.qza \
    --m-metadata-file metadata/metadata.tsv \
    --p-formula symptomatic \
    --p-reference-levels "symptomatic::no" \
    --p-conserve \
    --o-differentials ancombc/metabolites/symptomatic-diff.qza

qiime composition tabulate \
	--i-data ancombc/metabolites/symptomatic-diff.qza \
	--o-visualization ancombc/metabolites/symptomatic-diff.qzv
	
qiime composition da-barplot \
	--i-data ancombc/metabolites/symptomatic-diff.qza \
	--o-visualization ancombc/metabolites/symptomatic-da-barplot.qzv


### sample classification testing
mkdir sample_classification

# random forest
# infection status
qiime sample-classifier classify-samples-ncv \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--m-metadata-file metadata/metadata.tsv \
	--m-metadata-column infection_status \
	--p-cv 5 \
	--p-n-estimators 999 \
	--p-estimator 'RandomForestClassifier' \
	--p-parameter-tuning \
	--p-random-state 42 \
	--p-n-jobs 4 \
	--output-dir sample_classification/rf_infection_status

qiime sample-classifier confusion-matrix \
  --i-predictions sample_classification/rf_infection_status/predictions.qza \
  --i-probabilities sample_classification/rf_infection_status/probabilities.qza \
  --m-truth-file metadata/metadata.tsv \
  --m-truth-column infection_status \
  --o-visualization sample_classification/rf_infection_status/rf_infection_status_ncv_confusion_matrix.qzv

qiime metadata tabulate \
  --m-input-file sample_classification/rf_infection_status/feature_importance.qza \
  --o-visualization sample_classification/rf_infection_status/feature_importance.qzv
  
qiime sample-classifier heatmap \
  --i-table metabolites/full_metabolite_shared_scaled_table.qza \
  --i-importance sample_classification/rf_infection_status/feature_importance.qza \
  --m-sample-metadata-file metadata/metadata.tsv \
  --m-sample-metadata-column infection_status \
  --p-group-samples \
  --p-color-scheme 'viridis_r' \
  --p-feature-count 20 \
  --o-filtered-table sample_classification/rf_infection_status/important-feature-table-top-20.qza \
  --o-heatmap sample_classification/rf_infection_status/important-feature-heatmap-top-20.qzv
  
# symptomatic
qiime sample-classifier classify-samples-ncv \
	--i-table metabolites/full_metabolite_shared_scaled_table.qza \
	--m-metadata-file metadata/metadata.tsv \
	--m-metadata-column symptomatic \
	--p-cv 5 \
	--p-n-estimators 999 \
	--p-estimator 'RandomForestClassifier' \
	--p-parameter-tuning \
	--p-random-state 42 \
	--p-n-jobs 4 \
	--output-dir sample_classification/rf_symptomatic

qiime sample-classifier confusion-matrix \
  --i-predictions sample_classification/rf_symptomatic/predictions.qza \
  --i-probabilities sample_classification/rf_symptomatic/probabilities.qza \
  --m-truth-file metadata/metadata.tsv \
  --m-truth-column symptomatic \
  --o-visualization sample_classification/rf_symptomatic/rf_symptomatic_ncv_confusion_matrix.qzv

qiime metadata tabulate \
  --m-input-file sample_classification/rf_symptomatic/feature_importance.qza \
  --o-visualization sample_classification/rf_symptomatic/feature_importance.qzv
  
qiime sample-classifier heatmap \
  --i-table metabolites/full_metabolite_shared_scaled_table.qza \
  --i-importance sample_classification/rf_symptomatic/feature_importance.qza \
  --m-sample-metadata-file metadata/metadata.tsv \
  --m-sample-metadata-column symptomatic \
  --p-group-samples \
  --p-color-scheme 'viridis_r' \
  --p-feature-count 20 \
  --o-filtered-table sample_classification/rf_symptomatic/important-feature-table-top-20.qza \
  --o-heatmap sample_classification/rf_symptomatic/important-feature-heatmap-top-20.qzv



  
          