# Simulace driftu mrtvych tel delfinu v Jaderskem mori

# Autor: Martin Zirovnicky
# Uprava: integrace vstupnich bodu z CSV a vypocet metrik ze stredovych bodu bunek
# Aktualizace: pridani parametru "stav_tela" ovlivnujiciho delku driftu

install.packages(c("ncdf4", "ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "geosphere"))

library(ncdf4)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)

# 1. Nacteni dat o pevnine
pevnina <- ne_countries(scale = "medium", returnclass = "sf")

# 2. Nacteni dat o proudech

nc <- nc_open("proudy.nc")
uo <- ncvar_get(nc, "uo")
vo <- ncvar_get(nc, "vo")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

# 3. Nacteni vychozich bodu a vypocet centroidu bunek
sever <- read.csv("severni_mista.csv")
jih <- read.csv("jizni_mista.csv")

sever$lon <- (sever$left + sever$right) / 2
sever$lat <- (sever$top + sever$bottom) / 2
jih$lon <- (jih$left + jih$right) / 2
jih$lat <- (jih$top + jih$bottom) / 2

sever$skupina <- "sever"
jih$skupina <- "jih"
vychozi_body_raw <- rbind(sever, jih)

# 3a. Pretypovani na numericke hodnoty
vychozi_body_raw$lon <- as.numeric(vychozi_body_raw$lon)
vychozi_body_raw$lat <- as.numeric(vychozi_body_raw$lat)

# 3b. Pridani stavu tela pro kazdy bod (lze upravit dle realnych dat)
# Mozne hodnoty: "cerstve", "stredni", "pokrocile"
vychozi_body_raw$stav_tela <- sample(c("cerstve", "stredni", "pokrocile"), nrow(vychozi_body_raw), replace = TRUE)

# 4. Posunout body z pevniny
vychozi_body <- vychozi_body_raw
for (i in 1:nrow(vychozi_body)) {
  bod <- st_sfc(st_point(c(vychozi_body$lon[i], vychozi_body$lat[i])), crs = 4326)
  while (any(st_intersects(bod, pevnina, sparse = FALSE))) {
    vychozi_body$lat[i] <- vychozi_body$lat[i] - 0.05
    bod <- st_sfc(st_point(c(vychozi_body$lon[i], vychozi_body$lat[i])), crs = 4326)
  }
}

# 5. Funkce pro proud
najdi_validni_bod <- function(lon_pt, lat_pt, lon_grid, lat_grid, uo, vo) {
  lon_idx <- which.min(abs(lon_grid - lon_pt))
  lat_idx <- which.min(abs(lat_grid - lat_pt))
  u <- uo[lon_idx, lat_idx, 1]
  v <- vo[lon_idx, lat_idx, 1]
  return(list(u = u, v = v))
}

# 6. Simulace driftu
stav_tela_dny <- list(
  cerstve = 30,
  stredni = 15,
  pokrocile = 5
)

simuluj_drift_s_kontrolo <- function(lon0, lat0, lon, lat, uo, vo, pevnina, celkem_dni = 30, dt_hodiny = 6) {
  body <- data.frame(lon = numeric(), lat = numeric())
  aktualni_lon <- lon0
  aktualni_lat <- lat0
  body <- rbind(body, c(aktualni_lon, aktualni_lat))
  kroku <- celkem_dni * 24 / dt_hodiny
  dt <- dt_hodiny * 3600

  for (i in 1:kroku) {
    vysl <- najdi_validni_bod(aktualni_lon, aktualni_lat, lon, lat, uo, vo)
    u <- vysl$u
    v <- vysl$v
    if (is.na(u) || is.na(v)) break

    dx <- u * dt
    dy <- v * dt
    dx_deg <- dx / (111000 * cos(aktualni_lat * pi / 180))
    dy_deg <- dy / 111000
    aktualni_lon <- aktualni_lon + dx_deg
    aktualni_lat <- aktualni_lat + dy_deg

    bod <- st_sfc(st_point(c(aktualni_lon, aktualni_lat)), crs = 4326)
    if (any(st_intersects(bod, pevnina, sparse = FALSE))) break
    body <- rbind(body, c(aktualni_lon, aktualni_lat))
  }

  colnames(body) <- c("lon", "lat")
  return(body)
}

# 7. Simulace vsech trajektorii
trajektorie <- list()
for (i in 1:nrow(vychozi_body)) {
  stav <- vychozi_body$stav_tela[i]
  dny <- stav_tela_dny[[stav]]
  traj <- simuluj_drift_s_kontrolo(
    lon0 = vychozi_body$lon[i],
    lat0 = vychozi_body$lat[i],
    lon = lon, lat = lat, uo = uo, vo = vo, pevnina = pevnina,
    celkem_dni = dny
  )
  traj$bod <- i
  traj$skupina <- vychozi_body$skupina[i]
  trajektorie[[i]] <- traj
}

trajektorie_df <- do.call(rbind, trajektorie)
trajektorie_df$bod <- rep(1:length(trajektorie), sapply(trajektorie, nrow))
trajektorie_df$skupina <- rep(vychozi_body$skupina, sapply(trajektorie, nrow))

# 8. Vypocty metrik
trajektorie_df <- trajektorie_df[order(trajektorie_df$bod), ]
trajektorie_df$dist <- c(0, distHaversine(
  trajektorie_df[-nrow(trajektorie_df), c("lon", "lat")],
  trajektorie_df[-1, c("lon", "lat")]
))
trajektorie_df$dist[trajektorie_df$bod[-1] != trajektorie_df$bod[-nrow(trajektorie_df)]] <- 0
soucet_vzd <- aggregate(dist ~ bod + skupina, data = trajektorie_df, sum)

start_pos <- vychozi_body[, c("lon", "lat")]
final_pos <- aggregate(cbind(lon, lat) ~ bod, data = trajektorie_df, tail, 1)
displacement <- distHaversine(start_pos, final_pos[, c("lon", "lat")])

delka <- aggregate(bod ~ bod, data = trajektorie_df, length)
delka$hodin <- delka$bod * 6

metriky <- data.frame(
  bod = soucet_vzd$bod,
  skupina = soucet_vzd$skupina,
  total_distance_m = soucet_vzd$dist,
  direct_distance_m = displacement,
  drift_duration_h = delka$hodin
)

aggregate(. ~ skupina, data = metriky[, -1], FUN = mean)

# 9. Vykresleni trajektorii
ggplot() +
  geom_sf(data = pevnina, fill = "grey85", color = "black", size = 0.2) +
  geom_path(data = trajektorie_df, aes(x = lon, y = lat, group = bod, color = skupina), linewidth = 0.8) +
  geom_point(data = vychozi_body, aes(x = lon, y = lat, color = skupina), size = 2) +
  coord_sf(xlim = c(13, 20), ylim = c(40, 46), expand = FALSE) +
  labs(title = "Simulace driftu mrtvych tel delfinu v Jadranu", x = "Zemepisna delka", y = "Zemepisna sirka", color = "Skupina") +
  theme_minimal()
