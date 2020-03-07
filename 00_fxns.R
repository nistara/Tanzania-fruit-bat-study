# ==============================================================================
# *
# ==============================================================================
# add NAs when the length of the burst is below the numberSamples (792 default)
fix_acc_len = function(burst, numberSamples) {
    if(length(burst) < numberSamples) burst[ (length(burst) + 1):numberSamples ] = NA
    burst
}


# ==============================================================================
#  *
# ==============================================================================
# Function to assign flying/notflying activity using:
# 1. acc burst data, and
# 2. k-means

get_acc_activity = function(acc, numberSamples, axis) {

    cat("*****", acc$tagID[1], "*****\n")
    
    # create a matrix with the samples per acc burst
    accaxes = strsplit(acc$eobs.accelerations.raw, split = " ", fixed = TRUE) 
    accaxes = lapply(accaxes, fix_acc_len, numberSamples)
    
    accaxes2 = matrix(as.numeric(unlist(accaxes)), ncol = numberSamples, byrow = TRUE) 
    colnames(accaxes2) = paste0("V", 1:numberSamples)

    acc2 =data.frame(timestamp = acc$timestamp)

    xax = accaxes2[ , seq(1, numberSamples, 3)]
    yax = accaxes2[ , seq(2, numberSamples, 3)]
    zax = accaxes2[ , seq(3, numberSamples, 3)]

    acc2$var.x = matrixStats::rowVars(xax, na.rm = TRUE)
    acc2$var.y = matrixStats::rowVars(yax, na.rm = TRUE)
    acc2$var.z = matrixStats::rowVars(zax, na.rm = TRUE)

    # cluster the variances of the z and y axes via kmeans
    if(axis %in% "xy") {
    mykmeans = acc2 %>%
        select(var.x, var.y) %>%
        kmeans(2)
    }

    if(axis %in% "xz") {
        mykmeans = acc2 %>%
            select(var.x, var.z) %>%
            kmeans(2)
    }

    if(axis %in% "yz") {
        mykmeans = acc2 %>%
            select(var.y, var.z) %>%
            kmeans(2)
    }

    
    if(axis %in% "xyz") {
        mykmeans = acc2 %>%
            select(var.x, var.y, var.z) %>%
            kmeans(2)
    }

    acc2$cluster = mykmeans$cluster

    # determine which cluster is flying or not-flying
    agmn = aggregate(acc2$var.z, by = list(acc2$cluster), FUN = min) 
    fly_cluster = ifelse(agmn$x[1] > agmn$x[2], 1, 2)

    # add flying/not-flying to acc data
    acc2$activity = NA
    acc2$activity[ acc2$cluster %in% fly_cluster ] = "Flying"
    acc2$activity[ ! acc2$cluster %in% fly_cluster ] = "NotFlying"

    # return acc acivity + burst info
    list(acc = acc2,
         burst_info = list(xax = xax, yax = yax, zax = zax))

}


# ==============================================================================
# Get GPS behavior from axis classification
# ==============================================================================

get_acc_gps_behav = function(id, acc_activityL, gpsL) {
    cat("*****", id, "*****\n")
    acc = acc_activityL[[ id ]]$acc
    gps = gpsL[[ id ]]

    # get sunrise/sunset times and classify acc times
    gps_centroid = centroid(gps[ , c("location.long", "location.lat")])
    acc = get_day_night(acc, gps_centroid)
    
    # assign to each gps point:
    #   1. acc activity classification nearest in time
    #   2. time difference (sec) between the selected acc burst and gps point
    acc_near_gps = lapply(gps$timestamp, get_near_acc, acc) %>%
        data.table::rbindlist()

    gps[ , c("activity", "DayTime")]  =  acc[ acc_near_gps$minN,
                                             c("activity", "DayTime")]

    gps$tDiff = acc_near_gps$tDiff

    list(acc_behav = acc,
         gps_behav = gps,
         burst_info = acc_activityL[[ id ]]$burst_info)
    
}

# ==============================================================================
# *
# ==============================================================================

get_day_night = function(acc, gps_centroid) {
    
    sunrise = sunriset(gps_centroid, acc$timestamp,
                       proj4string = CRS("+proj=longlat +ellps=WGS84 	
                                    +datum=WGS84 +towgs84=0,0,0"), 
                       direction = c("sunrise"), POSIXct.out = TRUE)

    sunset = sunriset(gps_centroid, acc$timestamp,
                      proj4string = CRS("+proj=longlat +ellps=WGS84 	
                                    +datum=WGS84 +towgs84=0,0,0"), 
                      direction = c("sunset"), POSIXct.out = TRUE)

    acc$sunrise = sunrise$time
    acc$sunset = sunset$time

    acc$DayTime = ifelse(acc$timestamp > acc$sunrise & 
                         acc$timestamp < acc$sunset, "day", "night")

    acc

}



# ==============================================================================
# *
# ==============================================================================
get_near_acc = function(gps_time, acc) {
    tDiff = abs(acc$timestamp - gps_time)
    minN = which.min(tDiff)
    list(minN = minN, tDiff = as.difftime(tDiff[ minN ], units = "secs"))

}


# ==============================================================================
# * For alternating axis labels
# ==============================================================================
# ref: https://stackoverflow.com/a/42038321/5443003
label_fill <- function(orig, .offset=0, .mod=2, .fill=""){
    ## replace
    ii <- as.logical(
        ## offset==0 keeps first
        (1:length(orig)-1+.offset) %% .mod
    )
    orig[ii] <- .fill
    orig
}


# ==============================================================================
# * Break gps points by flying nights
# ==============================================================================

break_nights_fxn = function(df, write = FALSE) {
    df$dates = date(df$timestamp)
    date_seq = seq(min(df$dates), max(df$dates), by = "days")
    # days = as.Date(df$timestamp)

    if (length(unique(date_seq)) > 1) {

        time_cut = as.numeric(cut(df$timestamp,
                                  breaks = as.POSIXct(paste0(date_seq, " 17:00:00"),
                                                      tz = "Africa/Dar_es_Salaam")))

        # if(min(time_cut) > 1) {time_cut = time_cut - 1}
        df$time_cut = time_cut

        # Remove NAs, because this means gps fixes were taken before bat left
        #    roost, but didn't come back or tags fell off near the roost during
        #    the day
        df = df[ !is.na(df$time_cut), ]
    } else {
        df$time_cut = "<1"
    }
    df
}


# ==============================================================================
# make lines from GPS data
# ==============================================================================
# Helpful for visualization and mapping
make_lines = function(df) {
    lapply(split(df, df$time_cut), function(tc) {
        coordinates(tc) = c("location.long", "location.lat")
        pp = list(Lines(list(Line(coordinates(tc))), ID = "w"))
        SpatialLines(pp)
    })
}


# ==============================================================================
# Get cumulative distances
# ==============================================================================
get_dist_km = function(df, lon = "utm.easting", lat = "utm.northing") {
    coordinates(df) = c(lon, lat)
    pp = list(Lines(list(Line(coordinates(df))), ID = "w"))
    ps = SpatialLines(pp)
    # plot(ps)
    SpatialLinesLengths(ps, longlat = FALSE)/1000
}



# ==============================================================================
# # Map GPS points (foraging/day roosts, etc by foraging night
# ==============================================================================
map_gps = function(id, gps_list, gps_lines_list,
                   map_forage = FALSE,
                   map_roosts = FALSE,
                   map_notflying = FALSE) {

    cat("*** ", id, " ***\n")
    df = gps_list[[ id ]]
    df_lines_list = gps_lines_list[[ id ]]
    cols = c("#7F3C8D", "#11A579", "#3969AC", "#F2B701",
             "#E73F74", "#80BA5A", "#E68310", "#008695",
             "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99")
               
    map = leaflet() %>%
        addTiles()

    for(i in unique(df$time_cut)) {
        df_map = df[ df$time_cut %in% i, ]
        df_lines = df_lines_list[[ i ]]
        map = map %>%
            addCircles(data = df_map,
                       lng = df_map$location.long,
                       lat = df_map$location.lat,
                       radius = 3,
                       stroke = FALSE,
                       group = as.character(i),
                       popup = ~paste0("tagID: ", id, "<br>",
                                       "time: ", df_map$timestamp, "<br>",
                                       "time_cut: ", df_map$time_cut),
                       fillOpacity = ifelse(df$forage_pts, 1, 0.5),
                       opacity = ifelse(df$forage_time > 60, 1, 0.5),
                       color = cols[i]) %>%
            addPolylines(data = df_lines,
                         color = cols[i],
                         weight = 2,
                         group = as.character(i))
    }
    

    groups = as.character(unique(df$time_cut))
    
    if(map_roosts) {
        df_roosts = df[ df$roost %in% "new", ]

        if(nrow(df_roosts) > 0) {

            roost_colors = c("#66C5CC", "#F6CF71", "#F89C74",
                             "#DCB0F2", "#87C55F", "#9EB9F3")
            
            df_roosts$colors = cols[ df_roosts$time_cut ]
            df_roosts$colors = roost_colors[ as.factor(df_roosts$dates)]

            map = map %>% addCircleMarkers(data = df_roosts,
                                           lng = df_roosts$location.long,
                                           lat = df_roosts$location.lat,
                                           color = df$colors,
                                           radius = 5,
                                           popup = ~paste0("tagID: ",
                                                           id, "<br>",
                                                           "date: ",
                                                           df_roosts$roost_date, "<br>",
                                                           "time: ",
                                                           df_roosts$timestamp, "<br>",
                                                           "time_cut: ",
                                                           df_roosts$time_cut),
                                           group = "roosts")
            
        groups = c(groups, "roosts")

        }

    }

    if(map_notflying) {

        df_notflying = df[ df$activity %in% "NotFlying", ]
        
        map = map %>% addCircleMarkers(data = df_notflying,
                                       lng = df_notflying$location.long,
                                       lat = df_notflying$location.lat,
                                       color = "yellow",
                                       opacity = 1,
                                       fillOpacity = 1,
                                       radius = 3,
                                       stroke = FALSE,
                                       popup = ~paste0("tagID: ",
                                                       id, "<br>",
                                                       "time: ",
                                                       df_notflying$timestamp, "<br>",
                                                       "time_cut: ",
                                                       df_notflying$time_cut, "<br>",
                                                       "forage_time: ",
                                                       df$forage_time),
                                       group = "Not flying")

        groups = c(groups, "Not flying")
    }

    
    if(map_forage) {

        df_forage = df[ df$forage_pts, ]
        
        map = map %>% addCircleMarkers(data = df_forage,
                                       lng = df_forage$location.long,
                                       lat = df_forage$location.lat,
                                       color = "black",
                                       radius = 2,
                                       stroke = FALSE,
                                       popup = ~paste0("tagID: ",
                                                       id, "<br>",
                                                       "time: ",
                                                       df_forage$timestamp, "<br>",
                                                       "time_cut: ",
                                                       df_forage$time_cut, "<br>",
                                                       "forage_time: ",
                                                       df$forage_time),
                                       group = "forage")

        groups = c(groups, "forage")
    }

    map = map %>%
        addLayersControl(overlayGroups = groups,
                         options = layersControlOptions(collapsed = FALSE))

    map
}



# ------------------------------------------------------------------------------
# * Plot function
# ------------------------------------------------------------------------------

plot_gps = function(ID,
                    ID_group,
                    gps_forage,
                    tag_colors,
                    basemaps,
                    kilombero_wgs,
                    area_IDs,
                    p_height, p_width, ... ) {


    if(FALSE) {

        ID = "K5313"
        ID_group = "small"

        ID = "K5310"
        ID_group = "large"

        
    }
    

    ind_bat = gps_forage[[ ID ]]

    gp_IDs = area_IDs[[ ID_group ]]
    df_sf = gps_forage %>%
        rbindlist() %>%
        filter(tagID %in% gp_IDs) %>%
        st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)

    basemap = basemaps[[ ID_group ]]
    sc_info = extent(df_sf)
    
    col = "#5F4690"
    col = tag_colors[ ID ]
    
    new_roosts = which(ind_bat$roost %in% "new")

    ind_bat$time_cut = paste0("Night ", ind_bat$time_cut)

    g =

    ggmap(basemap) +
    geom_point(data = kilombero_wgs, aes(x = lon, y = lat),
                   color = "#08B8FF", size = 4, alpha = 0.65) +
    geom_path(data = ind_bat, aes(x = location.long, y = location.lat),
              color = col, size = 0.25) +
    geom_point(data = ind_bat, aes(x = location.long, y = location.lat),
               color = col, size = 0.5)
    # geom_point(data = kilombero_wgs, aes(x = lon, y = lat, shape = site),
               # color = "white", size = 0.5)

    if(length(new_roosts) > 0) {
        g = g + geom_point(data = ind_bat[ new_roosts[1], ],
                           aes(x = location.long, y = location.lat),
                           color = "#C33488", size = 4, alpha = 0.65)
    }

    g = g +
    geom_point(data = ind_bat[ ind_bat$forage_pts, ],
               aes(x = location.long, y = location.lat),
               color = "black", size = 1.5, alpha = 1) +
    geom_point(data = ind_bat[ ind_bat$forage_pts, ],
               aes(x = location.long, y = location.lat),
               color = "yellow", size = 1.25, alpha = 1) +
        scale_shape_manual(values = 8) +
        theme_void() +
        theme(
            plot.background = element_rect(fill = '#D1D6AB', colour = '#D1D6AB'),
            legend.position="none",
            plot.margin=grid::unit(c(1,1,1,1), "mm")
        ) +
        ggtitle(ID) 

    g_scale =  g +
        theme_ipsum(base_size = 10,
                    plot_title_size = 12,
                    grid = FALSE,
                    axis = FALSE, ticks = TRUE) +
        ggsn::scalebar(
                  x.min = sc_info[1] - 0.001,
                  x.max = sc_info[2] + 0.001,
                  y.min = sc_info[3] - ifelse(ID_group %in% "small", 0.01, 0.001),
                  y.max = sc_info[4] + 0.001,
                  dist = 4,
                  dist_unit = "km",
                  transform = TRUE,
                  model = "WGS84",
                  location = "bottomright",
                  height = ifelse(ID_group %in% "small", 0.02, 0.01),
                  st.dist = ifelse(ID_group %in% "small", 0.04, 0.025),
                  st.size = ifelse(ID_group %in% "small", 2.5, 2),
                  box.fill = c("white", "black"),
                  box.color = "black",
                  border.size = 0.25) +
        xlab("Longitude") +
        ylab("Latitude") +
        ggtitle(ID) 
    
    gf = g +
        facet_grid(. ~ time_cut) +
        theme_classic(base_size = 9) +
        theme(legend.position="none")

    fn = file.path(plot_out_dir, paste0(ID, ".png"))
    fn_scale = file.path(plot_out_dir, paste0(ID, "_scale.png"))
    fn_facet = file.path(plot_out_dir, "facet", paste0(ID, "_facet.png"))
    
    cowplot::save_plot(fn, g, base_height = 3, dpi = 600)
    cowplot::save_plot(fn_scale, g_scale, base_height = 5, dpi = 600)
    cowplot::save_plot(fn_facet, gf, base_height = 5, dpi = 600)
    
    print(g_scale)
    print(gf)
    cat("\n")
    
    g
    
}
