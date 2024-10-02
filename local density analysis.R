# =======================================================

# the distance_matrix function finds the difference from each point to every other point, 
# and returns a rectangular matrix (which means each distance is in the matrix twice, 
# once from point A to point B, and again from point B to point A)

distance_matrix <- function(location_data = locations, x = 1, y = 2){
  distance_m <- matrix(nrow = nrow(location_data), ncol = nrow(location_data))

  for (i in 1:nrow(location_data)) {
    for (j in 1:nrow(location_data)) {
      distance_m[i,j] <- sqrt((location_data[i,x] - location_data[j,x])^2 + 
                        (location_data[i,y] - location_data[j,y])^2)
    }
  }
  return(distance_m)
}

# =======================================================

# local_counts is a function to find the counts of points of different types within the 
# local neighborhood of each point -- points are included in the counts 
# of their own neighborhood. The output is a dataframe that includes the point location 
# and type of each point, as well as counts (by type) of points within the 
# specified radius of each point.

local_counts <- function(location_data = locations, radius){

  distance <- distance_matrix(location_data = location_data)
  output <- mutate(location_data, radius = radius) 
  
  type_list <- unique(location_data$type)
  type_list <- sort(type_list)
  
  # find the counts for artifacts of each type in the neighborhood of each point
  for (n in 1:length(type_list)) {
    current_type <- type_list[n]
  
    local_counts <- matrix(nrow = nrow(locations))
    for (i in 1:nrow(distance)) {
      in_distance <- 0
      for (j in 1:ncol(distance)) {
        if (locations$type[j] == current_type) {
          if (distance[i,j] <= radius) {
            in_distance <- in_distance + 1
          }  
        }
      }
      local_counts[i] <- in_distance
    }
    
 #this adds local counts for each type, one column at a time
    output <- cbind(output, local_counts)
    col_list <- colnames(output)
    colnames(output) <- c(col_list[1:ncol(output) - 1], paste0("count_", as.character(type_list[n])))

  }
# add the local counts for total points in the neighborhood
  local_counts <- matrix(nrow = nrow(location_data))
  
  for (i in 1:nrow(distance)) {
    in_distance <- 0
    for (j in 1:ncol(distance)) {
      if (distance[i,j] <= radius) {
        in_distance <- in_distance + 1
      }
    }
    local_counts[i] <- in_distance
  }
  output <- cbind(output, local_counts)
  
  col_list <- colnames(output)
  colnames(output) <- c(col_list[1:ncol(output) - 1], "count_total")
 
  return(output)
}


# =======================================================
# local_density is a function that calls the local_count function, 
# then calculates the area of the neighborhood defined by the specified 
# radius. The local_density is the count divided by the area of the 
# neighborhood. The function outputs a dataframe that includes the point 
# location and type of each point, as well as the density of points of 
# each type within the specified radius of each point.


local_density <- function(location_data = locations, radius) {
  counts <- local_counts(location_data, radius)
  local_densities <- counts
  type_list <- unique(local_densities$type)
  type_list <- sort(type_list)
  # reduce the same-type counts and the total counts for 
  # each point by 1 so that points don't count in calculating their
  # own local densities
  for (i in 1:nrow(local_densities)) {
    local_densities[i,ncol(local_densities)] <- local_densities[i,ncol(local_densities)] - 1
    for (j in 1:length(type_list))
    if (local_densities$type[i] == type_list[j]) {
      local_densities[i,j + 4] <- local_densities[i,j + 4] - 1
    }
  }
  area <- pi * radius^2
  col_list <- colnames(local_densities)
  for (i in 5:ncol(local_densities)) {
    col_list[i] <- substring(col_list[i],7)
    col_list[i] <- paste0("density_", col_list[i])
    for (j in 1:nrow(local_densities)) {
      local_densities[j,i] <- local_densities[j,i] / area
    }
  }
 colnames(local_densities) <- col_list 
  
 return(local_densities)  

}



# =======================================================
# glb_density calculates the global density of points of each type 
# (the number of points of the type divided by the total area)

glb_density <- function(location_data = locations, site_area) {
  type_list <- unique(location_data$type)
  type_list <- sort(type_list)
  global_density <- matrix(nrow = length(type_list) + 1)
  for (i in 1:length(type_list)) {
    temp <- filter(location_data, type == type_list[i])
    global_density[i] <- nrow(temp) / site_area
  }
  global_density[length(type_list) + 1] <- nrow(location_data) / site_area
  row.names(global_density) <- c(type_list, "total")
  return(global_density)
}


# =======================================================
#lda calculates the local density coefficient within and between 
# all the types in the dataset. The local density coefficient from Type A 
# to Type B is the mean density of points of Type B within the specified 
# radius of points of Type A, divided by the global density of points of 
# Type B (i.e., the number of points of Type B divided by the total area)

lda <- function(location_data = locations, radius, site_area) {
  for (n in 1:length(radius)) {
    densities <- local_density(location_data, radius[n])
    type_list <- unique(location_data$type)
    type_list <- sort(type_list)
    global_density <- glb_density(location_data, site_area)
    lda_matrix <- matrix(nrow = length(type_list) + 1, ncol = length(type_list) + 2)
    name_list <- c(as.character(type_list), "total")
    name_list_col <- c(paste0("ldc_", name_list), "radius")
    colnames(lda_matrix) <- name_list_col
    #need to be in a column rather than (or in addition to) 
    # as row names
    rownames(lda_matrix) <- name_list
  
      for (i in 1:(length(type_list) + 1)) {
      if (i <= length(type_list)){
        current_type <- filter(densities, type == type_list[i])
        for (j in 1:(length(type_list) + 1)) {
         lda_matrix[i,j] <- mean(current_type[,j + 4]) / global_density[j]
        }
      }
      else {
        for (j in 1:(length(type_list) + 1)) {
          lda_matrix[(length(type_list) + 1),j] <- mean(densities[,j + 4]) / global_density[j]
        }
      }
      
    }
    lda_matrix <- as.data.frame(round(lda_matrix, 2))
    type <- name_list
    lda_matrix <- cbind(type, lda_matrix)
    lda_matrix$radius <- radius[n]
    if (n == 1) {
      lda_out <- lda_matrix
    }
    else {
      lda_out <- rbind(lda_out, lda_matrix)
    }
  }

  return(lda_out)
}


# =======================================================

#code to test out the functions above with two example data files


#some functions use mutate() from the dplyr package
library(dplyr)


locations <- read.csv(file = file.choose())


# the data should be formatted as a data frame of artifact locations
# with x coordinates in the first column, y coordinates in the second, and artifact type
# in the third column

# this just ensures that the column names in the data file are consistent
colnames(locations) <- c("x", "y", "type")


# locations$type <- as.factor(locations$type) 

#site area for Kintigh's example file "LDEN.csv"
LDEN_site_area <- 154

# site area for "AZ_A1020_BLM_point_plots.csv"
A1020_site_area <- 2409

distance_test <- distance_matrix(location_data = locations)

counts_test <- local_counts(radius = 2)

local_density_test <- local_density(radius = 2)


global_density_test <- glb_density(site_area = A1020_site_area)


lda_test <- lda(locations, radius = cbind(1,2,5), site_area = AZ_A_10_20_area)
