library(GenomeInfoDbData)
data(speciesMap)

# crossmeta can get gene symbols for species in org_pkg or ens_spcs
load("~/Dropbox/GraduateSchool/PhD/LabWork/MetaOmic/transcriptomic/R/crossmeta/R/sysdata.rda")
# ens_spcs <- readRDS('/home/alex/Documents/Batcave/GEO/SRAdb/ens_spcs.rds')
ens_spcs <- list()

# already have taxon ids for species in ens_spcs
species <- setdiff(names(org_pkg), ens_spcs$scientific_name)

speciesMap <- speciesMap[speciesMap$species %in% species, ]
row.names(speciesMap) <- speciesMap$species

org_taxid <- c(speciesMap$taxon, speciesMap[species, 'taxon'])
names(org_taxid) <- c(speciesMap$species, species)

# devtools::use_data(gpl_bioc, homologene, sources, token, org_pkg, org_taxid, internal = TRUE, overwrite = TRUE)
