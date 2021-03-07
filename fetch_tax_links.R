args=commandArgs(trailingOnly = TRUE)
library(taxize)
blast <-read.delim(args[1], header = FALSE)

#Reading species list
print("Reading species list")
all_list <- blast[,5]
list <- as.character(unique(sort(all_list)))
sprintf("Processing %s species from %s hits", length(list), length(all_list))

#Check and correct the names via the Global Names Resolver (GNR) service provided by the Encyclopedia of Life.
print("Resolving names via GNR by EOL")
cor_list <- names(sapply(list, function(x) {
  res <-gnr_resolve(names=x, best_match_only = TRUE)
  res$matched_name
}))

#Fetching taxonomy information from ncbi.
print("Fetching taxonomy from ncbi")
taxonomy <- tax_name(query = cor_list, get = c("genus","family","order","class"), db = "ncbi")[,3:6]

#Fetching synonyms from worms.
print("Fetching synonyms from worms")
synonyms <- vector("character")
for (i in 1:length(cor_list)) {
  synonyms_i <- try(synonyms(cor_list[i], db ="worms", rows = 1))
  if (attr(synonyms_i,"class")=="try-error") {
    synonyms <- c(synonyms,"NA")
  } else if (!"scientificname" %in%names(synonyms_i[[1]])){
    synonyms <- c(synonyms, "NA")
  } else {
    synonyms <- c(synonyms, paste(synonyms_i[[1]]$scientificname,collapse = ";"))
  }
}

#Fetching Turkish common names from wikipages.
tur_names <- vector("character")
for (i in 1:length(cor_list)) {
  wikiinfo <- wikitaxa::wt_wikispecies(cor_list[i])
  if (length(wikiinfo)==0) {
    tur_names<-c(tur_names, "NA")
  } else if (length(wikiinfo$common_names$name[wikiinfo$common_names$language=="Türkçe"])==0) {
    tur_names<-c(tur_names, "NA")
  } else {
    tur_names<-c(tur_names, wikiinfo$common_names$name[wikiinfo$common_names$language=="Türkçe"])
  }
}

#Fetching WIKI links.
print("Fetching wiki links")
wiki <- vector("character")
for (i in 1:length(cor_list)) {
  wiki_i <- try(get_wiki(cor_list[i], rows = 1))
  if (attr(wiki_i,"class")=="try-error") {
    wiki <- c(wiki,"NA")
  } else if (attr(wiki_i,"match")=="not found") {
    wiki <- c(wiki,"NA")
  } else {
    wiki <- c(wiki, attr(wiki_i, "uri"))
  }
}

#Fetching BOLD links.
print("Fetching bold links")
bold <- vector("character")
for (i in 1:length(cor_list)) {
  bold_i <- try(get_boldid(cor_list[i], rows = 1))
  if (attr(bold_i,"class")=="try-error") {
    bold <- c(bold,"NA")
  } else if (attr(bold_i,"match")=="not found") {
    bold <- c(bold, "NA")
  } else {
    bold <- c(bold, attr(bold_i, "uri"))
  }
}

#Fetching WORMS links.
print("Fetching worms links")
worms <- vector("character")
for (i in 1:length(cor_list)) {
  worms_i <- try(get_wormsid(cor_list[i], rows = 1))
  if (attr(worms_i,"class")=="try-error") {
    worms <- c(worms,"NA")
  } else if (attr(worms_i, "match")=="not found") {
    worms <- c(worms,"NA")
  } else {
    worms <- c(worms, attr(worms_i, "uri"))
  }
}

#Fetching ITIS links.
print("Fetching itis links")
itis <- vector("character")
for (i in 1:length(cor_list)) {
  itis_i <- get_tsn(cor_list[i], rows = 1)
  if (attr(itis_i, "match")=="not found") {
    itis <- c(itis,"NA")
  } else {
    itis <- c(itis, attr(itis_i, "uri"))
  }
}

#Fetching COL links.
print("Fetching col links")
col <- vector("character")
for (i in 1:length(cor_list)) {
  col_i <- try(get_colid(cor_list[i], rows = 1))
  if (attr(col_i,"class")=="try-error") {
    col <- c(col,"NA")
  } else if (attr(col_i,"match")=="not found") {
    col <- c(col,"NA")
  } else {
    col <- c(col, attr(col_i, "uri"))
  }
}

#Fetching EOL links.
#print("Fetching eol links")
#eol <- vector("character")
#for (i in 1:length(cor_list)) {
#  eol_i <- eol_search(cor_list[i])
#  if (is.na(eol_i$name)) {
#    eol <- c(eol,"NA")
#  } else {
#    eol <- c(eol, eol_i$link[1])
#  }
#}

#Fetching GBIF links.
print("Fetching gbif links")
gbif <- vector("character")
for (i in 1:length(cor_list)) {
  gbif_i <- try(get_gbifid(cor_list[i], rows = 1))
  if (attr(gbif_i,"class")=="try-error") {
    gbif <- c(gbif, "NA")
  } else if (attr(gbif_i,"match")=="not found") {
    gbif <- c(gbif,"NA")
  } else {
    gbif <- c(gbif, attr(gbif_i,"uri"))
  }
}

#Fetching OTT links.
print("Fetching ott links")
ott <- vector("character")
for (i in 1:length(cor_list)) {
  ott_i <- try(get_tolid(cor_list[i], rows = 1))
  if (attr(ott_i,"class")=="try-error") {
    ott <- c(ott, "NA")
  } else if (attr(ott_i,"match")=="not found") {
    ott <- c(ott,"NA")
  } else {
    ott <- c(ott, attr(ott_i,"uri"))
  }
}

#Fetching IUCN links.
print("Fetching iucn links")
iucn <- vector("character")
for (i in 1:length(cor_list)) {
  iucn_i <- try(get_iucn(cor_list[i], key = "8ed2054663acc7e13d1f868697d0f3d99081af357df1c37a52e75abe5e3fe45d"))
  if (attr(iucn_i,"class")=="try-error") {
    iucn <- c(iucn, "NA")
  } else if (attr(iucn_i, "match")=="not found") {
    iucn <- c(iucn,"NA")
  } else {
    iucn <- c(iucn, attr(iucn_i, "uri"))
  }
}

mycolnames <- c("Code","BLAST Hit","Length","P. Ident","Scientific Name","WORMS Synonyms","Local Name","Genus","Family","Order","Class","WIKI Link","BOLD Link","WORMS Link","ITIS Link","COL Link","GBIF Link","OTT Link","IUCN Link")
tax_results <- cbind(synonyms, tur_names, taxonomy, wiki, bold, worms, itis, col, gbif, ott, iucn)
all_results <- vector("character")
for (i in 1:length(list)) {
  results_i <- cbind(blast[all_list %in% list[i],], tax_results[i,])
  all_results <- rbind(all_results, results_i)
}

write.table(all_results, file="results.xlsx", sep = "\t", quote = FALSE, row.names = FALSE, col.names = mycolnames)
