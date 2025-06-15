# Load necessary libraries
library(phyloseq)
library(dplyr)

# Load your initial R objects
load('R_objects/input.Rdata')

# 1. AGGLOMERATE & SELECT TOP TAXA
# Agglomerate taxa to the "Order" level
phy.phy <- tax_glom(phy, taxrank = "Order")

# Get the names of the top 50 taxa based on total counts
top_taxa_names <- names(sort(taxa_sums(phy.phy), decreasing = TRUE))[1:50]

# 2. TRANSFORM & SUBSET PHYLOSEQ OBJECT
# Calculate relative abundance
phy.rel <- transform_sample_counts(phy.phy, function(x) x / sum(x))

# Subset the relative abundance object to only the top taxa
phy.top <- prune_taxa(top_taxa_names, phy.rel)

# Assign the metadata to the phyloseq object

sample_data(phy.top) <- mtmpBETA

# 3. CREATE THE FINAL DATA FRAME
# Extract sample data and OTU table into clean data frames
sample_df <- data.frame(sample_data(phy.top), SampleID = sample_names(phy.top))
otu_df <- data.frame(otu_table(phy.top), SampleID = sample_names(phy.top))

# Join the two data frames by SampleID to create the wide-format table
datona <- left_join(sample_df, otu_df, by = "SampleID")

# 4. RENAME TAXON COLUMNS
# Create a lookup vector where names are Taxon IDs and values are Order names
# We use phyloseq:: to be explicit and avoid namespace conflicts
tax_lookup <- phyloseq::tax_table(phy.top)[, "Order"]
names(tax_lookup) <- phyloseq::taxa_names(phy.top)

# Get the current column names of the final data frame
current_names <- colnames(datona)

# Find which of the current names are taxon IDs that need to be replaced
names_to_replace <- intersect(current_names, names(tax_lookup))

# Replace the taxon IDs with the actual Order names
colnames(datona)[match(names_to_replace, current_names)] <- tax_lookup[names_to_replace]





# 5. SAVE THE FINAL OBJECT
# Assign the final cleaned data frame to your 'order' variable and save
order <- datona
save(order, file = "R_objects/metadataWithOrderTables.RData")