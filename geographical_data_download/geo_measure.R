# Load required libraries
library("countrycode")
library("CoordinateCleaner")
library("sf")
library("rgdal")
library("raster")
library("tidyverse")
library("rworldmap")
library("rgbif")
library("lwgeom")
library("exactextractr")

# Load your data
datos <- read_csv("total_species.csv", col_names = TRUE)
continent <- st_read("/home/fabianc.salgado/shared/Araneidae_phylogeny/geo_analysis/shapes/continent.shp")
join_islands <- st_read("/home/fabianc.salgado/shared/Araneidae_phylogeny/geo_analysis/shapes/islands_join.shp")
join_islands_sp <- as(join_islands, "Spatial")
continent_sp <- as(continent, "Spatial")
eco_reg <- st_read("/home/fabianc.salgado/shared/join_analysis_spiders/new_total/Ecoregions2017/Ecoregions2017.shp")
temperature <- raster("/home/fabianc.salgado/shared/join_analysis_spiders/new_total/Beck_KG_V1/Beck_KG_V1_present_0p0083.tif")
world <- st_as_sf(shapefile("/home/fabianc.salgado/shared/Araneidae_phylogeny/geo_analysis/World_Continents/World_Continents.shp"))
map_world <- getMap(resolution = "high")

# Load gbif data
gbif_data <- read_delim("gbif_total.csv", col_names = TRUE, delim = "\t")

# Create a list of species
list_species <- datos$species

# Create necessary folders
dir.create("results", showWarnings = FALSE)
dir.create("points", showWarnings = FALSE)

# Define function
geo_measure <- function(spe) {
  # Subset the GBIF data for the given species
  species_data <- gbif_data %>%
    filter(species == spe) %>%
    select(species, decimalLongitude, decimalLatitude, countryCode, gbifID, coordinateUncertaintyInMeters)

  # Drop rows with missing longitude or latitude
  species_data <- species_data %>%
    drop_na(decimalLongitude, decimalLatitude)

  # Check if there are enough reliable points
  if (nrow(species_data) == 0 || sum(species_data$coordinateUncertaintyInMeters) == 0) {
    stop("There are not trustable points on the GBIF for ", spe)
  }

  # Extract and clean the distribution data
  distri <- datos %>%
    filter(grepl(spe, species)) %>%
    pull(WSC_ISO3) %>%
    tolower() %>%
    strsplit(c("\\, | and | to |from |\\?|\\/|probably|possibly|introduced")) %>%
    lapply(sub, pattern = "\\(.*", replacement = "") %>%
    lapply(trimws)
  
  # Get ISO codes
  iso <- wscmap[wscmap[, 1] %in% unlist(distri), -1] %>%
    colSums() %>%
    names() %>%
    .[. > 0]
  
  # Remove specific regions
  excluded_regions <- c("SGS", "IOT")
  iso <- iso[!iso %in% excluded_regions]
  
  # Filter the world map to get relevant countries
  relevant_countries <- map_world[map_world$ISO3 %in% iso, ]
  
  # Calculate the area of relevant countries in square kilometers
  area_countries_km <- sum(area(relevant_countries)) / 1E6
  
  # Prepare spatial points
  species_points <- species_data %>% select(decimalLongitude, decimalLatitude) %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(relevant_countries))
  
  # Perform an overlay to identify points within relevant countries
  points_in_countries <- over(species_points, relevant_countries)

  # Create a new data frame with clean points
  cleaned_points <- species_data[points_in_countries$ScaleRank %in% 1:10, ]

  # Filter points with a resolution lower than 50 km
  cleaned_points <- cleaned_points %>%
    filter(coordinateUncertaintyInMeters / 1000 <= 50) %>%
    drop_na(coordinateUncertaintyInMeters)
  
  # Check if there are enough clean points
  if (nrow(cleaned_points) == 0) {
    stop("There are not enough points for ", spe, " to perform calculations")
  }

  # Save the clean points to a CSV file
  file_pattern <- file.path(getwd(), "points", paste(spe, ".csv", sep = ""))
  write_csv(cleaned_points, path = file_pattern, col_names = TRUE)
          # Calculate EcoRegions
        puntos_clean <- dat_cl %>%
          select(decimalLongitude, decimalLatitude) %>%
          drop_na(decimalLongitude) %>%
          drop_na(decimalLatitude) %>%
          mutate(
            decimalLatitude = as.numeric(decimalLatitude),
            decimalLongitude = as.numeric(decimalLongitude)
          )

        puntos_clean <- st_as_sf(puntos_clean, coords = c("decimalLongitude", "decimalLatitude"))
        puntos_clean <- st_set_crs(puntos_clean, st_crs(eco_reg))
        points_overlap <- over(as(puntos_clean, "Spatial"), as(eco_reg, "Spatial"))
        eco_reg_spe <- na.omit(unique(points_overlap$ECO_NAME))

        # Calculate Koppen climatic zones
        sp_temp <- raster::extract(temperature, as(puntos_clean, "Spatial"))
        climatic_zones <- length(unique(sp_temp[sp_temp > 0]))

        # Check if the points are on islands or continents
        puntos <- st_as_sf(dat_cl, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(join_islands))
        islands_points <- over(as(puntos, "Spatial"), as(join_islands, "Spatial"))
        puntos <- st_as_sf(dat_cl, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(continent))
        continent_points <- over(as(puntos, "Spatial"), as(continent, "Spatial"))

        # Categorize points as island, continent, or island_continent
        cat_island <- ifelse(
          length(which(is.na(islands_points$OBJECTID_1) == FALSE)) > 0,
          ifelse(
            length(which(is.na(continent_points$OBJECTID_1) == FALSE)) > 0,
            "island_continent",
            "island"
          ),
          ifelse(
            length(which(is.na(continent_points$OBJECTID_1) == FALSE)) > 0,
            "continent",
            NA
          )
        )

        # Calculate conditional WSC points for island and continent
        country_dist <- rgeos::gBuffer(country_dist, byid = TRUE, width = 0)
        join_islands_sp <- as(join_islands, "Spatial")
        crs(country_dist) <- crs(join_islands_sp)
        crs(continent_sp) <- crs(join_islands_sp)
        islands_wsc <- tryCatch(
          sp::over(country_dist, join_islands_sp),
          error = function(e) {
            options(sf_use_s2 = FALSE)
            sp::over(country_dist, join_islands_sp)
          }
        )

        continent_wsc <- tryCatch(
          sp::over(country_dist, continent_sp),
          error = function(e) {
            options(sf_use_s2 = FALSE)
            sp::over(country_dist, continent_sp)
          }
        )

        # Categorize based on WSC data for island and continent
        cat_island_wsc <- ifelse(
          length(which(is.na(islands_wsc$OBJECTID_1) == FALSE)) > 0,
          ifelse(
            length(which(is.na(continent_wsc$OBJECTID_1) == FALSE)) > 0,
            "island_continent",
            "island"
          ),
          ifelse(
            length(which(is.na(continent_wsc$OBJECTID_1) == FALSE)) > 0,
            "continent",
            NA
          )
        )

        # Calculate latitudinal range with points
        lat_range <- c(min(dat_cl$decimalLatitude), max(dat_cl$decimalLatitude))
        lat_range_total <- ifelse(
          lat_range[1] <= 0 && lat_range[2] >= 0,
          abs(lat_range[1]) + lat_range[2],
          abs(lat_range[2]) - abs(lat_range[1]
        ))

        # Calculate latitudinal range with the WSC info
        lat_range_wsc <- c(min(country_dist$LAT), max(country_dist$LAT))
        lat_range_wsc_total <- ifelse(
          lat_range_wsc[1] <= 0 && lat_range_wsc[2] >= 0,
          abs(lat_range_wsc[1]) + lat_range_wsc[2],
          abs(lat_range_wsc[2]) - abs(lat_range_wsc[1]
        ))

        # Calculate convex polygon
        sp_poly <- st_as_sf(dat_cl, coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84")
        sp_poly_df <- summarise(sp_poly, geometry = st_combine(geometry)) %>% st_convex_hull()
        sp_poly_df <- st_make_valid(sp_poly_df)

        # Intersect with world map
        exp_pol <- tryCatch(
          st_intersection(sp_poly_df, world[2]),
          error = function(e) {
            options(sf_use_s2 = FALSE)
            st_intersection(sp_poly_df, world[2])
          }
        )

        if (!exists("exp_pol") || typeof(exp_pol) != "list") {
          tmp_tabla <- tibble(
            species = spe,
            n_points = nrow(dat_cl),
            area_polygon = NA,
            area_countries_wsc = area_countries_km,
            lat_range = lat_range_total,
            lat_range_wsc = lat_range_wsc_total,
            cat_island = cat_island,
            cat_island_wsc = cat_island_wsc,
            centroid_lon = NA,
            centroid_lat = NA,
            eco_reg_points = length(eco_reg_spe),
            eco_reg_polygon = NA,
            temp_zones_points = climatic_zones,
            temp_zones_polygon = NA
          )
          archivo <- file.path(getwd(), "results", paste(spe, ".csv", sep = ""))
          write_csv(tmp_tabla, path = archivo, col_names = TRUE)
          cat("problems with", spe, "\n")
        } else {
          # Calculate the rest of your data and update the tibble
          # ...
          # Finally, save the data to the output CSV
          tmp_tabla <- tibble(
            species = spe,
            n_points = nrow(dat_cl),
            area_polygon = area_polygon,
            area_countries_wsc = area_countries_km,
            lat_range = lat_range_total,
            lat_range_wsc = lat_range_wsc_total,
            cat_island = cat_island,
            cat_island_wsc = cat_island_wsc,
            centroid_lon = unlist(centroide)[1],
            centroid_lat = unlist(centroide)[2],
            eco_reg_points = length(eco_reg_spe),
            eco_reg_polygon = length(eco_reg_country),
            temp_zones_points = climatic_zones,
            temp_zones_polygon = climatic_zones_wsc
          )
          archivo <- file.path(getwd(), "results", paste(spe, ".csv", sep = "")
          write_csv(tmp_tabla, path = archivo, col_names = TRUE)
          cat("all good with", spe, "\n")
        }
      }
    }
  }
}

# Run the function 

for(spe in datos$species){
    tryCatch(geo_measure(spe), error = function(e){})
}