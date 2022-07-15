#Load Libraries
library(raster)
library(mapview)
library(tidyverse)
library(dplyr)
library(sf)
library(spatstat)
library(sp)
library(maps)
library(mapdata)
library(ggpubr)
library(usmap)


###Map of USA###

#USA
usa = map_data("usa")

#Georgia
states = st_as_sf(map("state", plot = FALSE, fill = TRUE))
ga = states %>% filter(ID == "georgia")

#Eastern USA + Shade in Georgia + Square around ACC
usa_plot = ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), color = "black", fill = "white") + 
  geom_sf(data = ga, color = "black", fill = "grey") +
  geom_rect(aes(xmin = -84, xmax = -83.3, ymin = 33.5, ymax = 34), color = "black", fill = NA)  + 
  coord_sf(xlim = c(-94, -66), ylim = c(22, 50), expand = FALSE)  + 
  annotate(geom="text", x = -83.1, y = 32, label="Georgia", size = 2)  + 
  annotate(geom = "text", x = -83.6, y = 34.25, label="ACC", size = 2) +
  theme_bw() + 
  xlab("") + 
  ylab("")

tiff('./figures/usa_plot.tiff', units="in", width=4, height=5, res=600, compression = 'lzw')
plot(usa_plot)
dev.off()


  


#Map out Athens-Clarke County

ga = st_as_sf(map('county', plot=F, fill=T ))
ga = subset(ga, grepl('georgia', ga$ID))
acc = ga %>% separate(ID, c("state", "county"), sep =",") %>% filter(county == "clarke")

ggplot(data = acc) +
  geom_sf(color = "black")

#Layer onto a base map of Georgia 
states = st_as_sf(map("state", plot = FALSE, fill = TRUE))
ga = states %>% filter(ID == "georgia")



acc_plot = ggplot() +
  geom_sf(data = ga, fill = "white") + 
  geom_sf(data = acc, fill = "#10a53dFF") + 
  annotate(geom="text", x = -83, y = 32, label="Georgia", size = 5)  + 
  annotate(geom = "text", x = -84, y = 34.25, label="Athens-Clarke County") + 
  theme_void()



#Map out ACC Catchments
wrf_catchment = st_read("./data/wrf_catchment_simple/wrf_catchment_simple.shp")


wrf_catchment %>% ggplot() + 
  geom_sf(aes(fill = wrf), alpha=.80) + 
  scale_fill_manual(values=c("#ffcf20FF", "#2f9aa0FF", "#10a53dFF"), name="Catchment \n Region") +
  theme_minimal() 

#Now, make point observations for the plants
wrf_points = data.frame(wrf = c("A", "B", "C"), long = c(-83.360595, -83.390170, -83.318359), lat =  c(33.937635, 33.910025, 33.886025))
wrf_points_sf = st_as_sf(wrf_points, coords = c("long", "lat"), crs = 4269 )

wrf_labs = data.frame(wrf = c("A", "B", "C"), long = c(-83.360595, -83.390170, -83.318359), lat =  c(33.94763, 33.92002, 33.896025))
wrf_labs_sf = st_as_sf(wrf_labs, coords = c("long", "lat"), crs = 4269 )


wrf_plot = wrf_catchment %>% 
  ggplot() + 
  geom_sf(aes(fill = wrf), alpha = 0.8) + 
  scale_fill_manual(values=c("#440154FF", "#ffcf20FF", "#2f9aa0FF"), name="Catchment \n Region") +
  geom_sf_label(data = wrf_labs_sf, aes(label = paste("WRF", wrf, sep = "")), size = 5) +
  geom_sf(data = wrf_points_sf, size = 5, color = "black") +
  theme_bw() + 
  xlab("") + 
  ylab("")


tiff('./figures/wrf_plot.tiff', units="in", width=6, height=5, res=600, compression = 'lzw')
plot(wrf_plot)
dev.off()

