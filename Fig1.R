library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(scales)
library(viridis)

# 1. Read data
dat <- read_delim("mapdata2.csv", delim = NULL, show_col_types = FALSE)

# 2. Force longitude, latitude, Time to numeric, remove missing values
dat <- dat %>%
  mutate(
    Longitude_E = as.numeric(Longitude_E),
    Latitude_N  = as.numeric(Latitude_N),
    Time = as.numeric(Time)
  ) %>%
  drop_na(Longitude_E, Latitude_N, Type, Time) %>%
  filter(Longitude_E >= -180 & Longitude_E <= 180,
         Latitude_N >= -90 & Latitude_N <= 90)

# 3. Average time by longitude+latitude+vegetation type
point_dat <- dat %>%
  group_by(Longitude_E, Latitude_N, Type) %>%
  summarise(Time = mean(Time, na.rm = T), .groups = "drop")

# 4. Load world basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# 5. Plot: full map including Antarctica (Viridis color scheme)
ggplot() +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
            fill = "#e0e0e0", color = NA, inherit.aes = FALSE) +
  geom_sf(data = world, fill = "white", color = "gray50", linewidth = 0.3) +
  
  geom_point(
    data = point_dat,
    aes(x = Longitude_E, y = Latitude_N, fill = Time, shape = Type),
    size = 4.5, stroke = 0.7, alpha = 0.9
  ) +
  
  scale_shape_manual(
    values = c("Herbaceous" = 21, "Woody" = 23),
    name = "Plant Type"
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "Experiment Day",
    guide = guide_colorbar(title.position = "top", barwidth = 1.0, barheight = 8)
  ) +
  
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  
  scale_y_continuous(
    breaks = seq(-90, 90, by = 30),
    minor_breaks = seq(-90, 90, by = 15)
  ) +
  
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    minor_breaks = seq(-180, 180, by = 30)
  ) +
  
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    panel.grid = element_line(color = "#e0e0e0", linewidth = 0.2),
    
    plot.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
    
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
    legend.spacing = unit(2, "mm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_line(),
    # Fix: explicit ggplot2::margin + unit mm
    axis.title.y = element_text(margin = ggplot2::margin(r = 5, unit = "mm")),
    axis.title.x = element_text(margin = ggplot2::margin(t = 5, unit = "mm"))
  )

# Export vector graphic (for submission)
# ggsave("global_map.pdf", width = 16, height = 8, dpi = 600)