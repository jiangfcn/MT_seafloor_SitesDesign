function d = dis2geopoints(lon1,lat1,lon2,lat2)
    f=[6378.137,1/299.25722]; % eccentricity
    d = distance(lat1,lon1,lat2,lon2,f);
end
