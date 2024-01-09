# Work flow for obtain the geographical information for each species

load the necessary libraries to follow the workflow
```R
library("countrycode")
library("arakno")
library("rworldmap")
```

## Check that the names of the phylogeny are correct based on the World spider catalog

```R
# Get the world map in high resolution using the getMap function
map_world <- getMap(resolution = "high")

# Get the wsc data (assuming there's a function named wsc() that fetches data)
wsc()

# Check the scientific names of species in the species list
# Assuming species_theridiidae_information.csv contains a column named "species" with species names
species <- read_csv("species_theridiidae_information.csv", col_names = TRUE)
species$species <- gsub(pattern = "_", replacement = " ", species$species)

# Loop through each species name and use tryCatch to handle potential errors
for (spe in species$species) {
    tryCatch(
        checknames(spe),  # Assuming checknames is a function to validate scientific names
        error = function(e) {
            message(paste("Check the name of ", spe))
        }
    )
}


```

unfortenly the function has some errors, so you will need to change some names by hand

## World Spider Catalog Geographical information

Create a list with the geographical regions registered for each species in the world spider catalog. When there is a range this database (e.g. Colombia to Argentina), you have to replace this by the continent because the next step would be tricky

```R
# Read data from a CSV file named "data_filtered.csv" into a dataframe named 'species'
species <- read_csv("data_filtered.csv", col_names = TRUE)

# Create a new column named 'WSC_ISO3' and initialize it with NA values
species$WSC_ISO3 <- NA

# Replace underscores with spaces in the 'species' column
species$species <- gsub(pattern = "_", replacement = " ", species$species)

# Create a list with geographical regions from the World Spider Catalog (WSC)
for (spe in unique(species$species)) {
    # Extract distribution information for the current species from the WSC data (assumes wscdata is available)
    distri <- wscdata[wscdata$name == spe, 10]

    # Check if the distribution includes "Introduced"
    if (grepl("Introduced", distri) == TRUE) {
        distri <- gsub(pattern = "\\Introduced+.+", replacement = "", distri)  # Remove introduced locations
    }

    # Remove "Probably native to" from distribution
    distri <- gsub(pattern = "\\Probably native to ", replacement = "", distri)

    # Update the 'WSC_ISO3' column in the 'species' dataframe with the distribution information
    species[species$species == spe, "WSC_ISO3"] <- toString(distri)
}

# Write the modified 'species' dataframe to a new CSV file named "data_filtered_countries.csv"
write_delim(species, file = "data_filtered_countries.csv", delim = ",")
```

## Download GBIF data

```R
# Read data from a CSV file named "data_filtered_countries.csv" into a dataframe named 'datos'
datos <- read_csv("data_filtered_countries.csv", col_names = TRUE)

# Extract unique species names from 'datos' and get GBIF identifiers using taxize::get_gbifid_
list_species <- datos %>% pull(species) %>% taxize::get_gbifid_(method = "backbone")

# Create an empty tibble to store species names and corresponding GBIF keys
names_keys <- tibble(species = character(), key = double())

# Loop through each species in the 'list_species' list
for (spe in names(list_species)) {
    # Extract the GBIF key for each species with status "ACCEPTED"
    key <- list_species[[spe]] %>% filter(status == "ACCEPTED") %>% pull(usagekey)
    
    # Create a temporary tibble with species name and GBIF key
    tmp <- tibble(species = spe, key = key)
    
    # Combine the temporary tibble with the 'names_keys' tibble
    names_keys <- rbind(names_keys, tmp)
}

# GBIF credentials for downloading data
user <- "XXX"  # Your gbif.org username
pwd <- "XXX"    # Your gbif.org password
email <- "XXX@XXX"  # Your email

# The command that downloads occurrence data from GBIF using the obtained keys
# Note: To actually download this data, you need to go to the gbif website and login
occ_download(pred_in("taxonKey", names_keys$key), format = "SIMPLE_CSV", user = user, pwd = pwd, email = email)

```

Please note that the actual download must be done on the GBIF website after logging in.


## Run geo_measure.R to get the species data

__geo_measure.R__ performs geographical measurements and analyses on species data. For each species, the script generates CSV files containing the results of the geographical measurements. These files are saved in the "results" folder. The script will also generate CSV files with the "clean" geographical records stored in the "points" folder.

The input to run this files can be obtained here: 

datos: [here](https://github.com/fcsalgado/polymorphism_spider_rangesize/blob/main/geographical_data_download/data/data_filtered_countries.csv)
continent: [here](https://github.com/fcsalgado/polymorphism_spider_rangesize/blob/main/geographical_data_download/data/data_filtered_countries.csv)
join_islands: [here](https://github.com/fcsalgado/polymorphism_spider_rangesize/blob/main/geographical_data_download/data/data_filtered_countries.csv)
eco_reg: [here](https://ecoregions.appspot.com/)
temperature: [here](https://www.gloh2o.org/koppen/)
world: [here](https://github.com/fcsalgado/polymorphism_spider_rangesize/blob/main/geographical_data_download/data/data_filtered_countries.csv)
