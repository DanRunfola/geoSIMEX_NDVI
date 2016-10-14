# MacArther analysis using geoSIMEX
# ADM2-Based Analysis

if(run_analysis == TRUE)
{
library(maptools)
library(sp)
library(jsonlite)
library(doBy)
library(raster)
library(plyr)
library(stargazer)
source("/home/aiddata/Desktop/Github/geoSIMEX_NDVI/Analysis/Code/geoSIMEX.R")

#---------------------------------------------------#
##### * Settings * #####
#---------------------------------------------------#

source_dir <- "/home/aiddata/Desktop/Github/geoSIMEX_NDVI/Analysis"

data.subset.categories <- list(c("Cambodia", "Laos", "Myanmar", "Thailand", "Vietnam"))

data.subset.categories.i <- 3

# Parameters
analysis.year.begin <- 2001
analysis.year.end <- 2013
aid.year.begin <- 2005
aid.year.end <- 2010

pre.end <- aid.year.begin - 1
post.begin <- aid.year.end + 1

#---------------------------------------------------#
##### * Loading and Configuring Data * #####
#---------------------------------------------------#

##### Loading and Configuring MacArthur AidData #####
macArth.level_1a <- read.csv(paste(source_dir, "/Data/MacArthur_Geocoded_data/level_1a.csv", sep=""))

macArth.inf <- macArth.level_1a[macArth.level_1a$ad_sector_codes %in% c("210","220","230","320"),]
macArth.inf <- macArth.inf[macArth.inf$status %in% c("Completion","Implementation"),]
macArth.inf <- macArth.inf[(macArth.inf$transactions_end_year >= aid.year.begin) & 
                             (macArth.inf$transactions_end_year <= aid.year.end),]
macArth.inf$recipients_iso3 <- as.character(macArth.inf$recipients_iso3)
macArth.inf <- macArth.inf[macArth.inf$recipients_iso3 != "Unspecified",]
macArth.inf <- macArth.inf[!is.na(macArth.inf$project_id),]

# Aid Project Summary Statistics
macArth.inf <- macArth.inf[macArth.inf$recipients %in% c("Cambodia", "Laos", "Myanmar", "Thailand", "Viet Nam"),]

nrow(table(macArth.inf$project_id)) # Number of Projects
nrow(macArth.inf) # Number of Project Locations
sum(macArth.inf$even_split_commitments, na.rm=TRUE) # Total Commitments
table(macArth.inf$precision_code) # Precision Codes
table(macArth.inf$recipients)
# Years for aid projects are based on the transaction date of the project which extend from 2005 to 2010.
# The resulting dataset includes 65 discrete projects allocated across 222 project locations, totaling in $7,729,666,290 committed.
# Precision Codes
#  1  2  3  4  6  8 
# 90 12 63 34 18  5 

# Assinging PC6,8 Projects Lat/Long within Country (if dataset has PC6 and 8s)
for(iso3 in names(table(macArth.inf$recipients_iso3))){
  
  if(nrow(macArth.inf[macArth.inf$recipients_iso3 %in% iso3,][macArth.inf[macArth.inf$recipients_iso3 %in% iso3,]$precision_code %in% c(6, 8),]) > 0){
    macArth.inf[macArth.inf$recipients_iso3 %in% iso3,][macArth.inf[macArth.inf$recipients_iso3 %in% iso3,]$precision_code %in% c(6, 8),]$latitude <- macArth.inf[macArth.inf$recipients_iso3 %in% iso3,][macArth.inf[macArth.inf$recipients_iso3 %in% iso3,]$precision_code %in% c(1),]$latitude[1]
    macArth.inf[macArth.inf$recipients_iso3 %in% iso3,][macArth.inf[macArth.inf$recipients_iso3 %in% iso3,]$precision_code %in% c(6, 8),]$longitude <- macArth.inf[macArth.inf$recipients_iso3 %in% iso3,][macArth.inf[macArth.inf$recipients_iso3 %in% iso3,]$precision_code %in% c(1),]$longitude[1]
  } 
}

macArth.inf <- macArth.inf[!is.na(macArth.inf$latitude),]

# Defining Coordinates and Projecting
coordinates(macArth.inf) <- ~longitude+latitude
proj4string(macArth.inf) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

##### Extracting GADM Names to Aid Project Locations

# Load GADM Data
gadm.ecohotspot.adm2 <- readShapePoly(paste(source_dir, "/Data/ADM2/GADM_MacEcohotspotSubset_ADM2.shp", sep=""))

macArth.inf@data$NAME_2 <- as.character(extract(gadm.ecohotspot.adm2, macArth.inf)$NAME_2)
macArth.inf@data$NAME_1 <- as.character(extract(gadm.ecohotspot.adm2, macArth.inf)$NAME_1)
macArth.inf@data$NAME_0 <- as.character(extract(gadm.ecohotspot.adm2, macArth.inf)$NAME_0)

macArth.inf <- as.data.frame(macArth.inf@data)
# some ADM2 names are the same across different districts
macArth.inf$NAME_2 <- paste(macArth.inf$NAME_2, macArth.inf$NAME_1, macArth.inf$NAME_0)
macArth.inf$NAME_1 <- paste(macArth.inf$NAME_1, macArth.inf$NAME_0)

##### Loading and Configuring Covariate Data #####
dta.all.adm2 <- read.csv(paste(paste(source_dir, "/Data/extracts/covariates_allcountries_adm2.csv", sep="")))
dta.sea <- read.csv(paste(paste(source_dir, "/Data/extracts/sea.csv", sep="")))
dta.sea.ncc13 <- read.csv(paste(source_dir, "/Data/extracts/nighttimelights_2013/merge_sea_grid.csv", sep=""))
names(dta.sea.ncc13) <- c("ID", "ncc4_2013e")
dta.sea.merged <- merge(dta.sea, dta.sea.ncc13, by="ID")
dta.all.adm2 <- summaryBy(. ~ NAME_2 + NAME_1 + NAME_0, data=dta.sea.merged, keep.names=TRUE)
write.csv(dta.all.adm2, paste(source_dir, "/Data/extracts/covariates_SEA_adm2.csv", sep=""))

# some ADM2 names are the same across different districts; ensure are unique
dta.all.adm2$NAME_2 <- paste(dta.all.adm2$NAME_2, dta.all.adm2$NAME_1, dta.all.adm2$NAME_0)
dta.all.adm2$NAME_1 <- paste(dta.all.adm2$NAME_1, dta.all.adm2$NAME_0)

# Creating ADM IDs from ADM Names
dta.all.adm2$NAME_2.id <- as.numeric(as.factor(dta.all.adm2$NAME_2))
dta.all.adm2$NAME_1.id <- as.numeric(as.factor(dta.all.adm2$NAME_1))
dta.all.adm2$NAME_0.id <- as.numeric(as.factor(dta.all.adm2$NAME_0))
dta.all.adm2$NAME_2 <- as.character(dta.all.adm2$NAME_2)
dta.all.adm2$NAME_1 <- as.character(dta.all.adm2$NAME_1)
dta.all.adm2$NAME_0 <- as.character(dta.all.adm2$NAME_0)

##### Getting ADM Names and IDs into AidData Data

# Creating datasets with only adm names and merging that with aid project-location data
dta.all.adm2.namesIDs <- subset(dta.all.adm2, select=c(NAME_2, NAME_1, 
                                                       NAME_0, NAME_2.id, 
                                                       NAME_1.id, NAME_0.id))

# Merging
macArth.inf <- merge(macArth.inf, dta.all.adm2.namesIDs, by=c("NAME_2", "NAME_1", 
                                                              "NAME_0"), all.x=TRUE, all.y=FALSE)

# Some projects don't land in grid cell area 
macArth.inf <- macArth.inf[!is.na(macArth.inf$NAME_0.id),]

#---------------------------------------------------#
##### * Creating Variables * #####
#---------------------------------------------------#

dta.all.adm2$per_loss_2001_cuml <- dta.all.adm2$per_loss_2001
dta.all.adm2$per_loss_2002_cuml <- dta.all.adm2$per_loss_2001_cuml + dta.all.adm2$per_loss_2002
dta.all.adm2$per_loss_2003_cuml <- dta.all.adm2$per_loss_2002_cuml + dta.all.adm2$per_loss_2003
dta.all.adm2$per_loss_2004_cuml <- dta.all.adm2$per_loss_2003_cuml + dta.all.adm2$per_loss_2004
dta.all.adm2$per_loss_2005_cuml <- dta.all.adm2$per_loss_2004_cuml + dta.all.adm2$per_loss_2005
dta.all.adm2$per_loss_2006_cuml <- dta.all.adm2$per_loss_2005_cuml + dta.all.adm2$per_loss_2006
dta.all.adm2$per_loss_2007_cuml <- dta.all.adm2$per_loss_2006_cuml + dta.all.adm2$per_loss_2007
dta.all.adm2$per_loss_2008_cuml <- dta.all.adm2$per_loss_2007_cuml + dta.all.adm2$per_loss_2008
dta.all.adm2$per_loss_2009_cuml <- dta.all.adm2$per_loss_2008_cuml + dta.all.adm2$per_loss_2009
dta.all.adm2$per_loss_2010_cuml <- dta.all.adm2$per_loss_2009_cuml + dta.all.adm2$per_loss_2010
dta.all.adm2$per_loss_2011_cuml <- dta.all.adm2$per_loss_2010_cuml + dta.all.adm2$per_loss_2011
dta.all.adm2$per_loss_2012_cuml <- dta.all.adm2$per_loss_2011_cuml + dta.all.adm2$per_loss_2012
dta.all.adm2$per_loss_2013_cuml <- dta.all.adm2$per_loss_2012_cuml + dta.all.adm2$per_loss_2013
dta.all.adm2$per_loss_2014_cuml <- dta.all.adm2$per_loss_2013_cuml + dta.all.adm2$per_loss_2014

##### Pre-Trends #####
dta.all.adm2$airtemp.min.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$airtemp.max.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$airtemp.avg.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$precip.min.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$precip.max.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$precip.avg.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$lnyx_e.pre <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$lights.pre <- rep(0, nrow(dta.all.adm2))

for(i in analysis.year.begin:pre.end){
  dta.all.adm2$airtemp.min.pre <- dta.all.adm2$airtemp.min.pre + eval(parse(text = paste("dta.all.adm2$at41_",i,"m", sep="")))
  dta.all.adm2$airtemp.max.pre <- dta.all.adm2$airtemp.max.pre + eval(parse(text = paste("dta.all.adm2$at41_",i,"x", sep="")))
  dta.all.adm2$airtemp.avg.pre <- dta.all.adm2$airtemp.avg.pre + eval(parse(text = paste("dta.all.adm2$at41_",i,"e", sep="")))
  dta.all.adm2$precip.min.pre <- dta.all.adm2$precip.min.pre + eval(parse(text = paste("dta.all.adm2$pc41_",i,"m", sep="")))
  dta.all.adm2$precip.max.pre <- dta.all.adm2$precip.max.pre + eval(parse(text = paste("dta.all.adm2$pc41_",i,"x", sep="")))
  dta.all.adm2$precip.avg.pre <- dta.all.adm2$precip.avg.pre + eval(parse(text = paste("dta.all.adm2$pc41_",i,"e", sep="")))
  dta.all.adm2$lnyx_e.pre <- dta.all.adm2$lnyx_e.pre + eval(parse(text = paste("dta.all.adm2$lnyx_",i,"e", sep="")))
  dta.all.adm2$lights.pre <- dta.all.adm2$lights.pre + eval(parse(text = paste("dta.all.adm2$ncc4_",i,"e", sep="")))
}

dta.all.adm2$airtemp.min.pre <- dta.all.adm2$airtemp.min.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$airtemp.max.pre <- dta.all.adm2$airtemp.max.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$airtemp.avg.pre <- dta.all.adm2$airtemp.avg.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$precip.min.pre <- dta.all.adm2$precip.min.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$precip.max.pre <- dta.all.adm2$precip.max.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$precip.avg.pre <- dta.all.adm2$precip.avg.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$lnyx_e.pre <- dta.all.adm2$lnyx_e.pre / length(analysis.year.begin:pre.end)
dta.all.adm2$lights.pre <- dta.all.adm2$lights.pre / length(analysis.year.begin:pre.end)

##### Post-Trends #####
dta.all.adm2$airtemp.min.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$airtemp.max.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$airtemp.avg.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$precip.min.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$precip.max.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$precip.avg.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$lnyx_e.post <- rep(0, nrow(dta.all.adm2))
dta.all.adm2$lights.post <- rep(0, nrow(dta.all.adm2))

for(i in post.begin:analysis.year.end){
  dta.all.adm2$airtemp.min.post <- dta.all.adm2$airtemp.min.post + eval(parse(text = paste("dta.all.adm2$at41_",i,"m", sep="")))
  dta.all.adm2$airtemp.max.post <- dta.all.adm2$airtemp.max.post + eval(parse(text = paste("dta.all.adm2$at41_",i,"x", sep="")))
  dta.all.adm2$airtemp.avg.post <- dta.all.adm2$airtemp.avg.post + eval(parse(text = paste("dta.all.adm2$at41_",i,"e", sep="")))
  dta.all.adm2$precip.min.post <- dta.all.adm2$precip.min.post + eval(parse(text = paste("dta.all.adm2$pc41_",i,"m", sep="")))
  dta.all.adm2$precip.max.post <- dta.all.adm2$precip.max.post + eval(parse(text = paste("dta.all.adm2$pc41_",i,"x", sep="")))
  dta.all.adm2$precip.avg.post <- dta.all.adm2$precip.avg.post + eval(parse(text = paste("dta.all.adm2$pc41_",i,"e", sep="")))
  dta.all.adm2$lnyx_e.post <- dta.all.adm2$lnyx_e.post + eval(parse(text = paste("dta.all.adm2$lnyx_",i,"e", sep="")))
  dta.all.adm2$lights.post <- dta.all.adm2$lights.post + eval(parse(text = paste("dta.all.adm2$ncc4_",i,"e", sep="")))
}

dta.all.adm2$airtemp.min.post <- dta.all.adm2$airtemp.min.post / length(post.begin:analysis.year.end)
dta.all.adm2$airtemp.max.post <- dta.all.adm2$airtemp.max.post / length(post.begin:analysis.year.end)
dta.all.adm2$airtemp.avg.post <- dta.all.adm2$airtemp.avg.post / length(post.begin:analysis.year.end)
dta.all.adm2$precip.min.post <- dta.all.adm2$precip.min.post / length(post.begin:analysis.year.end)
dta.all.adm2$precip.max.post <- dta.all.adm2$precip.max.post / length(post.begin:analysis.year.end)
dta.all.adm2$precip.avg.post <- dta.all.adm2$precip.avg.post / length(post.begin:analysis.year.end)
dta.all.adm2$lnyx_e.post <- dta.all.adm2$lnyx_e.post / length(post.begin:analysis.year.end)
dta.all.adm2$lights.post <- dta.all.adm2$lights.post / length(post.begin:analysis.year.end)

##### Variables Change #####
dta.all.adm2$airtemp.min.change <- dta.all.adm2$airtemp.min.post - dta.all.adm2$airtemp.min.pre
dta.all.adm2$airtemp.max.change <- dta.all.adm2$airtemp.max.post - dta.all.adm2$airtemp.max.pre
dta.all.adm2$airtemp.avg.change <- dta.all.adm2$airtemp.avg.post - dta.all.adm2$airtemp.avg.pre
dta.all.adm2$precip.min.change <- dta.all.adm2$precip.min.post - dta.all.adm2$precip.min.pre
dta.all.adm2$precip.max.change <- dta.all.adm2$precip.max.post - dta.all.adm2$precip.max.pre
dta.all.adm2$precip.avg.change <- dta.all.adm2$precip.avg.post - dta.all.adm2$precip.avg.pre
dta.all.adm2$lnyx_e.change <- dta.all.adm2$lnyx_e.post - dta.all.adm2$lnyx_e.pre
dta.all.adm2$lights.change <- dta.all.adm2$lights.post - dta.all.adm2$lights.pre

# Change in Dependent Variable
dta.all.adm2$per_loss_cuml.change <- eval(parse(text = paste("dta.all.adm2$per_loss_",analysis.year.end,"_cuml", sep=""))) - eval(parse(text = paste("dta.all.adm2$per_loss_",pre.end,"_cuml", sep="")))

#---------------------------------------------------#
##### * Calculating Expected Value of Aid * #####
#---------------------------------------------------#

###### Transform Aid Commitments
# Log of Aid
#macArth.inf$even_split_commitments[is.na(macArth.inf$even_split_commitments)] <- 0
#macArth.inf$even_split_commitments <- log(macArth.inf$even_split_commitments+1)

# Number of Projects
macArth.inf$even_split_commitments <- 1

# Subset for low PC-only case; save all data.
macArth.inf.temp.lowPC <- macArth.inf[macArth.inf$precision_code %in% c(1,2,3),]


# Expected Values - (Using all precision codes)
dta.all.adm2$expected_aid <- expected_aid_ROI(aidData=macArth.inf, 
                                                   roiData=dta.all.adm2, 
                                                   probAidAssume=dta.all.adm2$Shape_Area, 
                                                   dollar_set=macArth.inf$even_split_commitments, 
                                                   aid.precision.code="precision_code", 
                                                   roi.pc1.name="NAME_2.id", 
                                                   roi.pc2.name="NAME_2.id", 
                                                   roi.pc3.name="NAME_2.id", 
                                                   roi.pc4.name="NAME_1.id", 
                                                   roi.pc5.name="NAME_1.id", 
                                                   roi.pc6.name="NAME_0.id", 
                                                   aid.pc1.centroid.name="NAME_2.id")



# Expected Values - (Using only precison codes 1, 2 and 3)
dta.all.adm2$expected_aid.lowPC <- expected_aid_ROI(aidData=macArth.inf.temp.lowPC, 
                                                   roiData=dta.all.adm2, 
                                                   probAidAssume=dta.all.adm2$Shape_Area, 
                                                   dollar_set=macArth.inf.temp.lowPC$even_split_commitments, 
                                                   aid.precision.code="precision_code", 
                                                   roi.pc1.name="NAME_2.id", 
                                                   roi.pc2.name="NAME_2.id", 
                                                   roi.pc3.name="NAME_2.id", 
                                                   roi.pc4.name="NAME_1.id", 
                                                   roi.pc5.name="NAME_1.id", 
                                                   roi.pc6.name="NAME_0.id", 
                                                   aid.pc1.centroid.name="NAME_2.id")



#---------------------------------------------------#
##### * Naive Models * #####
#---------------------------------------------------#

model.naive.expected_aid <- lm(per_loss_cuml.change ~ expected_aid + airtemp.min.change + airtemp.max.change + 
                                 airtemp.avg.change + precip.min.change + precip.max.change + precip.avg.change + 
                                 lnyx_e.change + lights.change, data=dta.all.adm2)

model.naive.expected_aid.lowPC <- lm(per_loss_cuml.change ~ expected_aid.lowPC + airtemp.min.change +
                                       airtemp.max.change + airtemp.avg.change + precip.min.change + 
                                       precip.max.change + precip.avg.change + lnyx_e.change + lights.change, 
                                     data=dta.all.adm2)

names(model.naive.expected_aid.lowPC$coefficients)[names(model.naive.expected_aid.lowPC$coefficients) == "expected_aid.lowPC"] <- "expected_aid"

#---------------------------------------------------#
##### * geoSIMEX Models * #####
#---------------------------------------------------#

model.geoSIMEX.expected_aid <- geoSIMEX(model=model.naive.expected_aid,
                                        geoSIMEXvariable="expected_aid",
                                        roiData=dta.all.adm2,
                                        aidData=macArth.inf,
                                        aid.amount="even_split_commitments",
                                        iterations=500,
                                        bins=4,
                                        roi.area="Shape_Area",
                                        roi.pc1.name="NAME_2.id", 
                                        roi.pc2.name="NAME_2.id", 
                                        roi.pc3.name="NAME_2.id", 
                                        roi.pc4.name="NAME_1.id", 
                                        roi.pc5.name="NAME_1.id", 
                                        roi.pc6.name="NAME_0.id", 
                                        aid.pc1.centroid.name="NAME_2.id",
                                        aid.precision.code="precision_code")

model.geoSIMEX.expected_aid$variance.model.imprecision$expected_aid[1] / sum(model.geoSIMEX.expected_aid$variance.model.imprecision$expected_aid)
model.geoSIMEX.expected_aid$variance.model.imprecision$expected_aid[2] / sum(model.geoSIMEX.expected_aid$variance.model.imprecision$expected_aid)


#---------------------------------------------------#
##### * Model Averaging Models * #####
#---------------------------------------------------#

model.modelAvg_randProb.expected_aid <- modelAverageRandProb(iterations=500, 
                                         aidData=macArth.inf, 
                                         roiData=dta.all.adm2, 
                                         aid.amount="even_split_commitments", 
                                         model=model.naive.expected_aid, 
                                         geoSIMEXvariable="expected_aid", binary=FALSE, 
                                         aid.precision.code="precision_code",
                                         roi.pc1.name="NAME_2.id", 
                                         roi.pc2.name="NAME_2.id", 
                                         roi.pc3.name="NAME_2.id", 
                                         roi.pc4.name="NAME_1.id", 
                                         roi.pc5.name="NAME_1.id", 
                                         roi.pc6.name="NAME_0.id", 
                                         aid.pc1.centroid.name="NAME_2.id")

#---------------------------------------------------#
##### * Stargazer Output * #####
#---------------------------------------------------#

# Expected Aid
CaseStudy_results_table <- stargazer(model.naive.expected_aid, model.naive.expected_aid.lowPC, model.modelAvg_randProb.expected_aid, model.geoSIMEX.expected_aid,
          title="Impact of Infrastructure Aid on Forest Loss", 
          label="seresults",
          covariate.labels=c("Aid",
                             "Air Temp (Max)",
                             "Air Temp (Min)",
                             "Air Temp (Avg)",
                             "Precip (Max)",
                             "Precip  (Min)",
                             "Precip  (Avg)",
                             "lnyx",
                             "Nighttime Lights"),
          omit.stat = c("f","ser","rsq","adj.rsq"),
          dep.var.caption = "Forest Loss",
          dep.var.labels   = c("",""),
          out="/home/aiddata/Desktop/Github/geoSIMEX_NDVI/Analysis/Results/CaseTable.tex",
          column.labels = c("Naive","Naive (Low PC)", "Model Avg", "geoSIMEX"),
          notes = "Standard errors in parentheses")

png("/home/aiddata/Desktop/Github/geoSIMEX_NDVI/images/geoSIMEX_plot.png")
plot(model.geoSIMEX.expected_aid, "expected_aid", name_variable="Aid")
dev.off()


}