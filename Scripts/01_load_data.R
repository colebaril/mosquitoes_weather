#' ---
#' title: "Mosquitoes - Load and Prep Data"
#' author: "Steffi LaZerte"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     theme: cosmo
#'     toc: true
#'     toc_float:
#'       collapsed: false
#' ---
#'
#' <!--- This style block is to tweak the formatting of the data tables in the report --->
#' <style>
#'   table {
#'     display: block;
#'     overflow: auto;
#'   }
#'
#'   blockquote {
#'     background-color: #88cabf;
#'  }
#'
#' </style>

#+ echo = FALSE
# RMarkdown -------------------------------------------------------------------

#' R Code to run and save this report:
#+ eval = FALSE
# rmarkdown::render(input = "Scripts/01_load_data.R",
#                   output_dir = "Results",
#                   output_file = paste0("01_load_data_", Sys.Date(), '.html'),
#                   envir = new.env())

#+ echo = FALSE
# Setup -------------------------------------------------------------------
#' # Setup

# renv::restore() # If using these scripts to restore dependencies (do not knit)

#+ setup, include = FALSE
knitr::opts_chunk$set(fig.width = 10, fig.asp = 0.5, out.width = "100%")

#' ## Load packages
#+ message = FALSE
library(tidyverse)
library(here)
library(lubridate)
library(assertr)

#' ## Create folders
dir.create(here("Data/Datasets"), showWarnings = FALSE)
dir.create(here("Results"), showWarnings = FALSE)

#' ## Load mosquitoes data

mosq <- read_csv(here("Data/Raw/Mosquito_Weather_Specific_Locations_Addedv3.csv"),
                 show_col_types = FALSE) %>%
  mutate(year = year(date)) %>%
  # Check data for formats, missing values and species
  assert(not_na, date, cdcweek, site, sitespecific, subsample, species) %>%
  assert(is.numeric, trapcount, ddm14, pt14, rha14, ttmax, ttmin, windspd) %>%
  assert(in_set(c("Ae. vexans", "Cq. perturbans", "Cx. tarsalis", "Oc. dorsalis")),
         species) %>%
  assert(in_set(seq(as_date("2020-05-01"), as_date("2021-09-01"), by = "1 day")), date)

#' ## Adding fog dates

fog <- data.frame(year = c(2020, 2021),
                  fog_date = as_date(c("2020-07-16", "2021-07-03")))

mosq <- mosq %>%
  mutate(year = year(date)) %>%
  left_join(fog, by = "year") %>%
  group_by(year) %>%
  mutate(after_fogged = site == "Brandon" & date >= fog_date) %>%
  ungroup()

#' Calculate lat lon for spatial autocorrelation

temp <- mosq %>%
  select(site) %>%
  mutate(site2 = if_else(site == "Portage", "Portage La Prairie",  site),
         site2 = paste0(site, ", MB, Canada")) %>%
  distinct()

geo <- geocode_OSM(temp$site2)
geo <- left_join(geo, temp, by = c("query" = "site2"))

mosq <- left_join(mosq, select(geo, "lat", "lon", "site"), by = "site")


#+ echo = FALSE
# Save Data ---------------------------------------------------------------
#' # Save Data
write_csv(mosq, here("Data/Datasets/mosquitoes.csv"))


#+ echo = FALSE
# Reproducible ---------------------------------------------------------------
#' # Reproducible

#' ## Data sets
DT::datatable(mosq, extensions = 'Buttons',
              options = list(dom = 'Bfrtip', buttons = c('csv', 'excel')))

#' ## Session Info
#+ R.options = list(width = 100)
devtools::session_info()

