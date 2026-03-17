@inline function iron_grazing(zoo::QualityDependantZooplankton, T, I, food_availability::NamedTuple, iron_availability::NamedTuple)
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensitivity
    p    = zoo.food_preferences
    food = keys(p)
    J    = zoo.specific_food_threshold_concentration
    K    = zoo.grazing_half_saturation
    food_threshold_concentration = zoo.food_threshold_concentration

    N = length(food)

    base_grazing_rate = g₀ * b ^ T

    total_food = sum(ntuple(n -> getproperty(food_availability, food[n]) * p[n], Val(N)))

    available_total_food = sum(ntuple(n -> max(zero(I), getproperty(food_availability, food[n]) - J) * p[n], Val(N)))

    concentration_limited_grazing = max(zero(I), available_total_food - min(available_total_food / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + total_food)

    total_specific_iron_grazing = sum(ntuple(n -> max(zero(I), getproperty(food_availability, food[n]) - J) * p[n] * getproperty(iron_availability, food[n]), Val(N))) * total_specific_grazing / (available_total_food + eps(zero(I)))

    return total_specific_iron_grazing * I
end

@inline function iron_grazing(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    food = prey_names(bgc, val_name)
    
    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    food_availability = extract_food_availability(bgc, i, j, k, fields, food)

    iron_ratios = extract_iron_availability(bgc, i, j, k, fields, food)

    return iron_grazing(zoo, T, I, food_availability, iron_ratios) 
end

@inline function iron_flux_feeding(zoo::QualityDependantZooplankton, T, I, sinking_iron_flux)
    g₀ = zoo.maximum_flux_feeding_rate
    b  = zoo.temperature_sensitivity

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_iron_flux 

    return total_specific_flux_feeding * I
end

@inline function iron_flux_feeding(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = edible_iron_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    return iron_flux_feeding(zoo, T, I, sinking_flux)
end
