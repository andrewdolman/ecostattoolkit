# Import lake test dataset

library(tidyverse)

readxl::excel_sheets("data-raw/Datensatz_Seen_Jahrgänge bis 2014_TN_TP_EQR.xls")

lake_mp_1 <- readxl::read_excel("data-raw/Datensatz_Seen_Jahrgänge bis 2014_TN_TP_EQR.xls", sheet = 2) 

lake_mp <- lake_mp_1 %>% 
  mutate(
    Unique_ID = OWK,
    EQR = `MP-EQR`,
    P = `TP_µg/L`,
    Exclude_P = 0,
    N = `TN_µg/L`,
    Exclude_N = 0,
    Group = `WRRL-Typ`) %>% 
  select(Unique_ID, EQR, P, Exclude_P, N, Exclude_N, Group) %>% 
  group_by(Unique_ID) %>% 
  summarise_all(mean) %>% 
  mutate(BioClass = cut(EQR,
                        breaks = seq(0, 1, by = 0.2),
                        labels = rev(c("High", "Good", "Moderate", "Poor", "Bad"))))



lake_pp_1 <- readxl::read_excel("data-raw/Datensatz_Seen_Jahrgänge bis 2014_TN_TP_EQR.xls", sheet = 3) 

lake_pp <- lake_pp_1 %>% 
  separate(col=MessstellennameJahr, into=c("Messstellenname", "Jahr"), sep = -5) %>% 
    mutate(
    Unique_ID = Messstellenname,
    EQR = EQR_PSI,
    P = TPSais,
    Exclude_P = 0,
    N = TNSais_µg,
    Exclude_N = 0,
    Group = `Typ-Gruppe gleiche Reftrophie`) %>% 
  select(Unique_ID, EQR, P, Exclude_P, N, Exclude_N, Group) %>% 
  group_by(Unique_ID) %>% 
  summarise_all(mean) %>% 
  mutate(BioClass = cut(EQR,
                        breaks = seq(0, 1, by = 0.2),
                        labels = rev(c("High", "Good", "Moderate", "Poor", "Bad"))))


devtools::use_data(lake_mp, lake_pp)

write.csv(lake_mp, paste0("inst/extdata/lake-test-data/", "lake_mp", ".csv"), row.names = FALSE, fileEncoding = "Latin1")
write.csv(lake_pp, paste0("inst/extdata/lake-test-data/", "lake_pp", ".csv"), row.names = FALSE, fileEncoding = "Latin1")


