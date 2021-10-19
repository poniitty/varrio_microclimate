##########################################################################
# THIS CODE MAKES SOME QUALITY CHECKS IF THE SITE AND TOMST IDs LOOKS FINE
#

library(tidyverse)

# List binary and command files to be removed from repository if also data file exists

f <- c(list.files("data", pattern = "binary_", recursive = T, full.names = T),
       list.files("data", pattern = "command_", recursive = T, full.names = T))

for(i in f){ if(file.exists(gsub("binary_","data_",i)) | file.exists(gsub("command_","data_",i))){
  unlink(i)
} else {
  print(paste0("DATA FILE MISSING!!! ", i))
} 
}
# If no printed messages then no problems


###########################################################################
# Check Tomst ID-numbers from last year data
maxdt <- read_csv("data/reading_times_2020.csv") %>% 
  mutate(site = as.character(site))

f <- list.files("data", pattern = "data_", recursive = T, full.names = T)

fi <- data.frame(file = f)

fi$site <- toupper(unlist(lapply(fi$file, function(x) strsplit(x, "/")[[1]][3])))

fi <- fi[order(fi$site),]

fi$tomst_id <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",strsplit(x, "/")[[1]][4]), "_")[[1]][1])))

fi %>% group_by(tomst_id) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(tomst_id) -> doubled_ids
fi %>% filter(tomst_id %in% doubled_ids) # check for weird things!!! Good if none

# Check if more than one data file in a folder
fi %>% group_by(site) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(site) -> doubled_sites
fi %>% filter(site %in% doubled_sites) # check for weird things!!! Good if none

###########################################################################
# Update the file list

f <- list.files("data", pattern = "data_", recursive = T, full.names = T)

fi <- data.frame(file = f)

fi$site <- toupper(unlist(lapply(fi$file, function(x) strsplit(x, "/")[[1]][3])))

fi <- fi[order(fi$site),]

fi$tomst_id <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",strsplit(x, "/")[[1]][4]), "_")[[1]][1])))

fi %>% group_by(tomst_id) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(tomst_id) -> doubled_ids
fi %>% filter(tomst_id %in% doubled_ids) # check for weird things!!! Good if none

fi %>% group_by(site) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(site) -> doubled_sites
fi %>% filter(site %in% doubled_sites) # check for weird things!!! Good if none
# Looks good now!

#######################################################################
# Check if missing sites in 2021 data
all <- full_join(fi, maxdt %>% rename(tomst_id_20 = tomst_id))

# Check for duplicate sites
all %>% group_by(site) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(site) -> doubled_sites
all %>% filter(site %in% doubled_sites) # No, Good!

# Non-matching sites
all %>% filter(!complete.cases(.))

# No sites that occur only in 2021 data

#sites 13, 19, 39, 55 in 2020 data but not in 2021
all %>% filter(tomst_id == 94194329) # No such Tomst ids in 2021 so this is fine!
all %>% filter(tomst_id == 94194328) # No such Tomst ids in 2021 so this is fine!
all %>% filter(tomst_id == 94194299) # No such Tomst ids in 2021 so this is fine!
all %>% filter(tomst_id == 94194325) # No such Tomst ids in 2021 so this is fine!

# For sites 13, 19, 39, 55 find 2020 data and copy to repository
f2 <- list.files("C:/Users/OMISTAJA/OneDrive - University of Helsinki/R_Projects/microclim_suomi/raw_field_data",
                 pattern = "data_", recursive = T, full.names = T)

# Copy site 13 data from last year data
f2[grepl("/varrio/tomst_varrio_2020/13", f2)]
dir.create(paste0(getwd(), "/data/tomst/13"))
file.copy(f2[grepl("/varrio/tomst_varrio_2020/13", f2)],
          paste0(getwd(), "/data/tomst/13/data_94194329_0.csv"))

# Copy site 19 data from last year data
f2[grepl("/varrio/tomst_varrio_2020/19", f2)]
dir.create(paste0(getwd(), "/data/tomst/19"))
file.copy(f2[grepl("/varrio/tomst_varrio_2020/19", f2)],
          paste0(getwd(), "/data/tomst/19/data_94194328_0.csv"))

# Copy site 39 data from last year data
f2[grepl("/varrio/tomst_varrio_2020/39", f2)]
dir.create(paste0(getwd(), "/data/tomst/39"))
file.copy(f2[grepl("/varrio/tomst_varrio_2020/39", f2)],
          paste0(getwd(), "/data/tomst/39/data_94194299_0.csv"))

# Copy site 55 data from last year data
f2[grepl("/varrio/tomst_varrio_2020/55", f2)]
dir.create(paste0(getwd(), "/data/tomst/55"))
file.copy(f2[grepl("/varrio/tomst_varrio_2020/55", f2)],
          paste0(getwd(), "/data/tomst/55/data_94194325_0.csv"))


########################################################################################
# Update file list

f <- list.files("data", pattern = "data_", recursive = T, full.names = T)

fi <- data.frame(file = f)

fi$site <- toupper(unlist(lapply(fi$file, function(x) strsplit(x, "/")[[1]][3])))

fi <- fi[order(fi$site),]

fi$tomst_id <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",strsplit(x, "/")[[1]][4]), "_")[[1]][1])))

fi %>% group_by(tomst_id) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(tomst_id) -> doubled_ids
fi %>% filter(tomst_id %in% doubled_ids) # check for weird things!!! Good if none

fi %>% group_by(site) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(site) -> doubled_sites
fi %>% filter(site %in% doubled_sites) # check for weird things!!! Good if none
# Looks good still!!!

#######################################################################
# Check if Tomst ids match between years
all <- full_join(fi, maxdt %>% rename(tomst_id_20 = tomst_id))

# Check for duplicate sites
all %>% group_by(site) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(site) -> doubled_sites
all %>% filter(site %in% doubled_sites) # No, Good!

all %>% filter(tomst_id == tomst_id_20)
all %>% filter(tomst_id != tomst_id_20)
# All seems to match nicely!!!!!!!!!!


# Good to go and read the data


