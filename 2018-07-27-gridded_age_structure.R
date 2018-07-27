# Grid interpolation of European regional age structure
# Jonas Schöley
# 2018-07-27

# Init --------------------------------------------------------------------

library(eurostat)      # eurostat data
library(rnaturalearth) # worldwide map data
library(tidyverse)     # tidy data transformation
library(lubridate)     # date and time support
library(sf)            # simple features GIS

# Download data -----------------------------------------------------------

# download the data on pop counts by age at NUTS-3 level
euro_age <-
  get_eurostat('demo_r_pjanaggr3', stringsAsFactors = FALSE) %>%
  # filter NUTS-3, 2015, total population
  filter(sex == 'T', str_length(geo) == 5,  year(time) == 2015,
         age %in% c('Y_LT15', 'Y15-64', 'Y_GE65', 'TOTAL')) %>%
  select(age, geo, values) %>%
  spread(age, values) %>%
  rename(total = TOTAL, age65plus = Y_GE65, age0to15 = Y_LT15, age15to65 = `Y15-64`)

# download geospatial data for NUTS-3 regions
# and project to crs 3035
euro_nuts3_sf <-
  get_eurostat_geospatial(output_class = 'sf',
                          resolution = '60', nuts_level = 3) %>%
  st_transform(crs = 3035)

# download geospatial data for European, Asian and African countries
eura <-
  ne_countries(continent = c('europe', 'asia', 'africa'), returnclass = 'sp') %>%
  fortify() %>%
  rename(longitude = long, latitude = lat)

# divide the european continent into a 100 by 100 km grid
euro_grid <-
  st_make_grid(euro_nuts3_sf, crs = 3035, what = 'polygons', cellsize = 100*1e3)

# SANITY CHECK: are the grid cell areas identical
range(st_area(euro_grid))

# interpolate population counts by age over the grid
# and extract grid centroids
euro_grid_age_sf <-
  left_join(euro_nuts3_sf, euro_age, by = 'geo') %>%
  gather(key, value, total, age65plus, age0to15, age15to65) %>%
  group_by(key) %>%
  do(
  # area weighted grid interpolation preserving totals
    select(., value) %>%
      st_interpolate_aw(to = euro_grid, extensive = TRUE) %>%
      st_centroid()
  ) %>%
  # projection flag got lost, reset it
  st_as_sf(crs = 3035) %>%
  spread(key, value)

# SANITY CHECK: do the totals make sense?
# should be approx 0
euro_grid_age_sf %>%
  # calculate population shares by grid cell
  mutate(diff = age0to15+age15to65+age65plus - total) %>%
  {range(.$diff)}

# average European age structure in 2015
euro_center <- with(na.omit(euro_age),
                    c(p_age0to15 = sum(age0to15)/sum(total),
                      p_age15to65 = sum(age15to65)/sum(total),
                      p_age65plus = sum(age65plus)/sum(total)))

# add population shares by age and
# differences from European average
euro_grid_age_sf <-
  euro_grid_age_sf %>%
  mutate(p_age0to15 = age0to15/total,
         p_age15to65 = age15to65/total,
         p_age65plus = age65plus/total,
         d_age0to15 = p_age0to15-euro_center['p_age0to15'],
         d_age15to65 = p_age15to65-euro_center['p_age15to65'],
         d_age65plus = p_age65plus-euro_center['p_age65plus'])

# prepare clean data frame
euro_grid_age_df <-
  euro_grid_age_sf %>%
  # project back to spherical (lon-lat) coordinates
  st_transform(4326) %>%
  # add lon-lat coordinates of grid centroids
  cbind(st_coordinates(.)) %>%
  # convert to standard data frame
  as_data_frame() %>%
  select(-geometry) %>%
  rename(grid_id = Group.1, longitude = X, latitude = Y) %>%
  arrange(desc(total))

# SANITY CHECK: do the differences make sense?
# should be approx 0
euro_grid_age_df %>%
  mutate(sum = d_age0to15 + d_age15to65 + d_age65plus) %>%
  {range(.$sum)}

write_csv(euro_grid_age_df, 'euro_grid_age.csv')

# Gridded population age structure ----------------------------------------

# color code the compositional deviation from european avg
euro_grid_age_df$col <-
  tricolore::Tricolore(euro_grid_age_df, legend = FALSE,
                       p1 = 'p_age65plus', p2 = 'p_age15to65', p3 = 'p_age0to15',
                       center = rev(prop.table(euro_center)), spread = 2.9,
                       contrast = .5, lightness = 1, chroma = 1, hue = 2/12)

# bubble-grid-map
euro_grid_age_df %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_polygon(aes(group = group), data = eura, color = 'white', fill = 'grey95') +
  geom_point(aes(col = col, size = total), show.legend = FALSE) +
  coord_map(xlim = c(-30, 50), ylim = c(35, 70), projection = 'lambert', 43, 62) +
  scale_size_area(
    'Population size',
    max_size = 8,
    breaks = c(1e3, 1e4, 1e5, 5e5),
    labels = c('1,000', '10,000', '100,000', '500,000')) +
  scale_color_identity() +
  theme_void() +
  theme(legend.position = c(0.83, 0.7)) +
  labs(caption = 'Data: Eurostat')

