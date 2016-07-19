library(parallel)
result_list <- mclapply(1, function(n)source("/home/aiddata/Desktop/Github/geoSIMEX_NDVI/Analysis/Code/MacArthur_Analysis.R"))

for(i in 1:length(result_list)){
  if(i == 1){
    results.all.df <- result_list[[1]]$value
  } else {
    results.all.df <- rbind(results.all.df, result_list[[1]]$value)
  }
}

# Generating Time to Save Results
op <- options(digits.secs = 3)
randN <- Sys.time()
randN <- gsub("-", "", randN, fixed = TRUE)
randN <- gsub(" ", "", randN, fixed = TRUE)
randN <- gsub(":", "", randN, fixed = TRUE)
randN <- gsub(".", "", randN, fixed = TRUE)

save(results, file=paste("/home/aiddata/Desktop/Github/geoSIMEX_NDVI/Analysis/Results/McArthur_MC_",randN,".Rda",sep=""))

