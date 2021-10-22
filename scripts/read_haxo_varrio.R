# Read HAXO logger data

library(tidyverse)
library(lubridate)
library(data.table)
library(cowplot)
library(zoo)

# List Haxo files
f <- list.files("data/haxo/", pattern = ".csv$", full.names = T)

fi <- data.frame(file = f,
                 site = parse_number(unlist(lapply(f, function(x) strsplit(x, "/")[[1]][3]))))

# read old HAXO data
old <- read_csv("data/haxo_data_corrected_2020.csv")

old %>% 
  mutate(datetime = with_tz(datetime, tzone = "Etc/GMT-2")) %>% 
  group_by(site) %>% 
  summarise(min_time = max(datetime)) -> old_times

old_times <- full_join(old_times, read_csv("data/reading_times_2021.csv") %>% 
                         mutate(maxdt = with_tz(maxdt, tzone = "Etc/GMT-2")) %>% 
                         rename(max_time = maxdt))

old_times %>% filter(!complete.cases(.))
old_times[old_times$site == 42, "max_time"] <- as_datetime("2021-07-09 14:00:00", tz = "Etc/GMT-2")

fi %>% filter(!site %in% unique(old$site)) # Good no new site id's

fi %>% left_join(., old_times) -> fi

fi <- fi[order(fi$site),]

# Loop through all the files and combine data
df <- data.frame()
for(i in fi$file){
  #i <- "./haxo_varrio_2020/54.csv"
  # i <- fi$file[3]
  print(i)
  
  d <- fread(i, header = T)
  
  if(NCOL(d) > 6){
    names(d) <- paste0("V",1:ncol(d))
    d %>% mutate(at = as.numeric(paste(V4,V5,sep = ".")),
                 arh = as.numeric(paste(V6,V7,sep = "."))) %>% 
      rename(Date = V2,
             Time = V3,
             Type = V8) %>% 
      select("V1","Date","Time","at","arh","Type") -> d
  } else {
    names(d) <- c("V1","Date","Time","at","arh","Type")
  }
  
  d$site <- fi[which(fi$file == i),"site"]
  
  # Prepare the Haxo data
  d %>% select(-V1,-Type) %>%
    mutate(datetime = ymd_hms(paste(Date, Time), tz = "Etc/GMT-3")) %>% 
    mutate(datetime = round_date(datetime, "hour")) %>% 
    mutate(datetime = with_tz(datetime, tzone = "Etc/GMT-2")) %>% 
    select(-Date,-Time) %>% 
    mutate(at = as.numeric(gsub(",",".",at))) %>% 
    mutate(arh = as.numeric(gsub(",",".",arh))) %>% 
    relocate(site, datetime) -> d
  
  if(d %>% slice(1) %>% pull(datetime) %>% hour() %% 2L == 0){
    d %>% mutate(datetime = datetime + hours(1)) -> d
  }
  
  # Filter bad data
  d %>% filter(datetime >= fi[which(fi$file == i),"min_time"],
               datetime <= fi[which(fi$file == i),"max_time"]) -> d
  
  df <- rbind.data.frame(df, d)
}

df %>% arrange(site, datetime) -> df

df %>% filter(datetime > "2019-07-01",
              datetime < "2021-08-01") -> df

sites <- unique(df$site)

pdf("visuals/Haxo_graphs.pdf", 12, 10)
for(i in sites){
  #i <- sites[3]
  print(i)
  
  df %>% filter(site == i) %>%
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "at"), col = "cornflowerblue") +
    theme_minimal() +
    ylab("Temperature") + xlab("Date")+
    scale_y_continuous(limits = c(-30, 35))+
    ggtitle(i) -> GG1
  
  df %>% filter(site == i) %>%
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "arh"), col = "cornflowerblue") +
    theme_minimal() +
    ylab("Air humidity") + xlab("Date")+
    scale_y_continuous(limits = c(0, 100)) -> GG2
  
  
  print(plot_grid(plotlist = list(GG1,GG2), nrow = 2))
  
}
dev.off()


####################################################################
# SPOT AND MARK ERRORS AND SUSPICIOUS ONES
# Site 51 reported to have fallen to the ground

df %>% group_by(datetime) %>% 
  summarise(med_at = median(at, na.rm = T),
            med_arh = median(arh, na.rm = T)) %>% 
  ungroup() %>% as.data.table() -> df_median

df %>% left_join(., df_median) -> df

df2 <- data.frame()
for(i in sites){
  
  print(i)
  
  df %>% filter(site == i) %>% 
    mutate(cor85_at = rollapply(., 85 ,function(x) cor(as.numeric(x[,"at"]),
                                                       as.numeric(x[,"med_at"])),
                                by.column=FALSE, fill = NA, partial = T)) %>% 
    mutate(cor85_arh = rollapply(., 85 ,function(x) cor(as.numeric(x[,"arh"]),
                                                        as.numeric(x[,"med_arh"])),
                                 by.column=FALSE, fill = NA, partial = T)) %>% 
    mutate(date = as_date(datetime)) %>% 
    mutate(rollsd_at = rollapply(at, width=61, FUN=sd, fill = NA, partial = T)) %>% 
    mutate(rollsd_arh = rollapply(arh, width=61, FUN=sd, fill = NA, partial = T)) %>% 
    mutate(abser_at = abs(at-med_at),
           abser_arh = abs(arh-med_arh)) %>% 
    mutate(rollabser_at = rollapply(abser_at, width=61, FUN=mean, fill = NA, partial = T)) %>% 
    mutate(rollabser_arh = rollapply(abser_arh, width=61, FUN=mean, fill = NA, partial = T)) -> temp
  
  
  
  df2 <- bind_rows(df2, temp)
  
}

df2 %>% group_by(site, date) %>% 
  summarise(mean_at = mean(at),
            min_at = min(at),
            max_at = max(at),
            sd_at = sd(at),
            sd_arh = sd(arh),
            meansd_at = mean(rollsd_at),
            meansd_arh = mean(rollsd_arh),
            mean_arh = mean(arh),
            min_arh = min(arh),
            cor85_at = min(cor85_at),
            cor85_arh = min(cor85_arh),
            rollabser_at = mean(rollabser_at),
            rollabser_arh = mean(rollabser_arh)) %>% 
  ungroup() %>% as.data.table() -> daily

daily %>% mutate(haxo_probl = 0) -> daily

###################################################
# Look into problematic sites one by one
# 51
df %>% filter(site == 51) %>% 
  filter(datetime > "2020-09-01",
         datetime < "2020-11-01") %>% 
  ggplot(aes_string(x="datetime")) +
  geom_line(aes_string(y = "at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("Temperature") + xlab("Date")
df %>% filter(site == 51) %>% 
  filter(datetime > "2020-09-01",
         datetime < "2020-11-01") %>% 
  ggplot(aes_string(x="datetime")) +
  geom_line(aes_string(y = "arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("Temperature") + xlab("Date")

daily %>% 
  mutate(haxo_probl = ifelse(site == 51 & date > "2020-09-20", 1, haxo_probl)) -> daily


###############################################################
# WHICH ARE UNDER SNOW

daily %>% 
  mutate(cor85_at = ifelse(is.na(cor85_at), 0, cor85_at),
         cor85_arh = ifelse(is.na(cor85_arh), 0, cor85_arh)) %>% 
  mutate(cor85_at = ifelse(!is.finite(cor85_at), 0, cor85_at),
         cor85_arh = ifelse(!is.finite(cor85_at), 0, cor85_arh)) -> daily

daily %>% filter(site == 51) -> temp

temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "max_at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("max_at") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "min_at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("min_at") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "sd_at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("sd_at") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "sd_arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("sd_arh") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "meansd_at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("meansd_at") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "meansd_arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("meansd_arh") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "mean_arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("mean_arh") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "min_arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("min_arh") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "cor85_at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("cor85_at") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "cor85_arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("cor85_arh") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "rollabser_at"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("rollabser_at") + xlab("Date") + ggtitle(unique(temp$site))
temp %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "rollabser_arh"), col = "cornflowerblue") +
  theme_minimal() +
  ylab("rollabser_arh") + xlab("Date") + ggtitle(unique(temp$site))

count_consecutives <- function(x){
  
  x <- ifelse(is.na(x), 666, x)
  
  x1 <- rle(x)
  
  x2 <- rep(x1$lengths, times = x1$lengths)
  
  return(x2)
}

daily %>% mutate(snow = 0) %>%
  mutate(snow = ifelse(min_arh > 90 & 
                         mean_arh > 92 &
                         meansd_arh < 1.5 &
                         meansd_at < 1.8 & 
                         sd_arh < 1.5 &
                         sd_at < 1.5 &
                         min_at > -15 &
                         max_at < 0.5 &
                         cor85_at < 0.8 &
                         cor85_arh < 0.8 &
                         rollabser_at > 1 &
                         rollabser_arh > 2, 2, snow)) %>% 
  mutate(haxo_probl = ifelse(haxo_probl == 1, haxo_probl, snow)) -> daily

daily %>% group_by(site) %>% 
  mutate(cons_count = count_consecutives(haxo_probl)) %>% 
  ungroup() %>% 
  mutate(haxo_probl = ifelse(haxo_probl == 0 & cons_count < 4, 2, haxo_probl)) %>% 
  group_by(site) %>% 
  mutate(cons_count = count_consecutives(haxo_probl)) %>% 
  ungroup() %>% 
  mutate(haxo_probl = ifelse(haxo_probl == 2 & cons_count < 7, 0, haxo_probl)) %>% 
  group_by(site) %>% 
  mutate(lag_probl = rollapply(haxo_probl, width=6, FUN=max, fill = NA, partial = T, align = "left")) %>% 
  mutate(lead_probl = rollapply(haxo_probl, width=3, FUN=max, fill = NA, partial = T, align = "right")) %>% 
  mutate(haxo_probl = ifelse(haxo_probl == 0 & lag_probl == 2, 2, haxo_probl)) %>% 
  mutate(haxo_probl = ifelse(haxo_probl == 0 & lead_probl == 2, 2, haxo_probl)) %>% 
  ungroup() -> daily2

daily2 %>% filter(site == 51) %>% 
  ggplot(aes_string(x="date")) +
  geom_line(aes_string(y = "mean_at"), col = "cornflowerblue") +
  geom_point(aes(x = date, y = mean_at),
             data = daily2 %>% filter(site == 51) %>% filter(haxo_probl == 2),
             col = "red") +
  geom_point(aes(x = date, y = mean_at),
             data = daily2 %>% filter(site == 51) %>% filter(haxo_probl == 1),
             col = "black") +
  theme_minimal() +
  ylab("min_arh") + xlab("Date")

df %>% mutate(date = as_date(datetime)) %>% 
  left_join(., daily2 %>% select(site, date, haxo_probl)) -> df2

pdf("visuals/Haxo_graphs_marked.pdf", 12, 10)
for(i in sites){
  #i <- sites[4]
  print(i)
  
  df2 %>% filter(site == i) %>% 
    mutate(at1 = as.numeric(ifelse(haxo_probl == 1, at, NA))) %>% 
    mutate(at2 = as.numeric(ifelse(haxo_probl == 2, at, NA))) %>% 
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "at"), col = "black") +
    geom_line(aes_string(y = "at2"), col = "blue") +
    geom_line(aes_string(y = "at1"), col = "red") +
    theme_minimal() +
    ylab("min_arh") + xlab("Date") +
    scale_y_continuous(limits = c(-30, 35))+
    ggtitle(i) -> GG1
  
  df2 %>% filter(site == i) %>% 
    mutate(arh1 = as.numeric(ifelse(haxo_probl == 1, arh, NA))) %>% 
    mutate(arh2 = as.numeric(ifelse(haxo_probl == 2, arh, NA))) %>% 
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "arh"), col = "black") +
    geom_line(aes_string(y = "arh2"), col = "blue") +
    geom_line(aes_string(y = "arh1"), col = "red") +
    theme_minimal() +
    ylab("Air humidity") + xlab("Date")+
    scale_y_continuous(limits = c(0, 100)) -> GG2
  
  
  print(plot_grid(plotlist = list(GG1,GG2), nrow = 2))
  
}
dev.off()

# MANUAL CORRECTION

df2 %>% mutate(haxo_probl = ifelse(site == 2, 0, haxo_probl)) -> df2
df2 %>% mutate(haxo_probl = ifelse(site == 47 & date > "2021-02-01" & date < "2021-04-01", 2, haxo_probl)) -> df2

fwrite(df2 %>% select(-med_at, -med_arh), "output/haxo_data.csv")
