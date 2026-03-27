module Iron

export SimpleIron, IronTendencyArgs

using Oceananigans.Units
using Oceananigans.Grids: znode, Center

using OceanBioME.Models.PISCESModel: PISCES, anoxia_factor

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: 
    aggregation_of_colloidal_iron, degradation

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    iron_scavenging, iron_scavenging_rate, bacterial_iron_uptake,
    specific_degradation_rate, edible_flux_rate, edible_iron_flux_rate

using OceanBioME.Models.PISCESModel.Phytoplankton: uptake

using OceanBioME.Models.PISCESModel.Zooplankton: 
    non_assimilated_iron, upper_trophic_dissolved_iron,
    bacteria_concentration, bacteria_activity

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: free_iron

include("simple_iron.jl")

@inline function ligand_concentration(iron::SimpleIron, DOC)
    Lₜᵐᵃˣ = iron.maximum_ligand_concentration
    Lₜ = iron.dissolved_ligand_ratio * DOC - Lₜᵐᵃˣ

    return max(Lₜᵐᵃˣ, Lₜ)
end

@inline function ligand_aggregation(iron::SimpleIron, Fe, DOC, T, scavenging_rate)
    total_ligand_concentration = ligand_concentration(iron, DOC)
    free_iron_concentration = free_iron(iron, Fe, DOC, T)
    excess_iron = max(zero(Fe), Fe - total_ligand_concentration)

    return iron.excess_scavenging_enhancement * scavenging_rate * excess_iron * free_iron_concentration
end

@inline function iron_tendency(iron::SimpleIron, args::IronTendencyArgs)
    λFe = iron_scavenging_rate(args.minimum_iron_scavenging_rate,
                               args.load_specific_iron_scavenging_rate,
                               args.POC,
                               args.GOC,
                               args.CaCO₃,
                               args.PSi)
    Fe′ = free_iron(iron, args.Fe, args.DOC, args.T)

    food_availability = (; P = args.P,
                           D = args.D,
                           POC = args.POC,
                           Z = args.Z)
    iron_availability = (; P = args.PFe / (args.P + eps(zero(args.P))),
                           D = args.DFe / (args.D + eps(zero(args.D))),
                           POC = args.SFe / (args.POC + eps(zero(args.POC))),
                           Z = args.microzooplankton_iron_ratio)

    Bact = bacteria_concentration(args.microzooplankton_bacteria_concentration,
                                  args.mesozooplankton_bacteria_concentration,
                                  args.maximum_bacteria_concentration,
                                  args.bacteria_concentration_depth_exponent,
                                  args.z,
                                  args.zₘₓₗ,
                                  args.zₑᵤ,
                                  args.Z,
                                  args.M)
    LBact = bacteria_activity(args.doc_half_saturation_for_bacterial_activity,
                              args.nitrate_half_saturation_for_bacterial_activity,
                              args.ammonia_half_saturation_for_bacterial_activity,
                              args.phosphate_half_saturation_for_bacterial_activity,
                              args.iron_half_saturation_for_bacterial_activity,
                              args.NH₄,
                              args.NO₃,
                              args.PO₄,
                              args.Fe,
                              args.DOC)

    small_particle_iron_remineralisation = degradation(Val(:SFe),
                                                       specific_degradation_rate(args.base_breakdown_rate,
                                                                                 args.particle_temperature_sensitivity,
                                                                                 args.O₂,
                                                                                 args.T,
                                                                                 args.first_anoxia_threshold,
                                                                                 args.second_anoxia_threshold),
                                                       args.SFe)
    grazing_waste = non_assimilated_iron(args.microzooplankton_iron_ratio,
                                         args.microzooplankton_non_assimilated_fraction,
                                         args.microzooplankton_maximum_grazing_rate,
                                         args.microzooplankton_temperature_sensitivity,
                                         args.microzooplankton_preference_for_p,
                                         args.microzooplankton_preference_for_d,
                                         args.microzooplankton_preference_for_z,
                                         args.microzooplankton_preference_for_poc,
                                         args.microzooplankton_specific_food_threshold_concentration,
                                         args.microzooplankton_grazing_half_saturation,
                                         args.microzooplankton_food_threshold_concentration,
                                         args.microzooplankton_minimum_growth_efficiency,
                                         args.microzooplankton_maximum_flux_feeding_rate,
                                         args.mesozooplankton_iron_ratio,
                                         args.mesozooplankton_non_assimilated_fraction,
                                         args.mesozooplankton_maximum_grazing_rate,
                                         args.mesozooplankton_temperature_sensitivity,
                                         args.mesozooplankton_preference_for_p,
                                         args.mesozooplankton_preference_for_d,
                                         args.mesozooplankton_preference_for_z,
                                         args.mesozooplankton_preference_for_poc,
                                         args.mesozooplankton_specific_food_threshold_concentration,
                                         args.mesozooplankton_grazing_half_saturation,
                                         args.mesozooplankton_food_threshold_concentration,
                                         args.mesozooplankton_minimum_growth_efficiency,
                                         args.mesozooplankton_maximum_flux_feeding_rate,
                                         args.T,
                                         args.Z,
                                         args.M,
                                         food_availability,
                                         iron_availability,
                                         args.sinking_flux,
                                         args.sinking_iron_flux)
    upper_trophic_waste = upper_trophic_dissolved_iron(args.mesozooplankton_minimum_growth_efficiency,
                                                       args.mesozooplankton_non_assimilated_fraction,
                                                       args.mesozooplankton_iron_ratio,
                                                       args.mesozooplankton_temperature_sensitivity,
                                                       args.mesozooplankton_quadratic_mortality,
                                                       args.T,
                                                       args.M)
    phytoplankton_iron_uptake = uptake(args.nano_exudated_fraction,
                                       args.nano_maximum_iron_ratio,
                                       args.nano_half_saturation_for_iron_uptake,
                                       args.nano_threshold_for_size_dependency,
                                       args.nano_size_ratio,
                                       args.nano_minimum_ammonium_half_saturation,
                                       args.nano_minimum_nitrate_half_saturation,
                                       args.nano_optimal_iron_quota,
                                       args.nano_base_growth_rate,
                                       args.nano_temperature_sensitivity,
                                       args.diatom_exudated_fraction,
                                       args.diatom_maximum_iron_ratio,
                                       args.diatom_half_saturation_for_iron_uptake,
                                       args.diatom_threshold_for_size_dependency,
                                       args.diatom_size_ratio,
                                       args.diatom_minimum_ammonium_half_saturation,
                                       args.diatom_minimum_nitrate_half_saturation,
                                       args.diatom_optimal_iron_quota,
                                       args.diatom_base_growth_rate,
                                       args.diatom_temperature_sensitivity,
                                       Val(:Fe),
                                       args.T,
                                       args.Fe,
                                       args.NO₃,
                                       args.NH₄,
                                       args.PO₄,
                                       args.Si,
                                       args.Si′,
                                       args.P,
                                       args.PChl,
                                       args.PFe,
                                       args.D,
                                       args.DChl,
                                       args.DFe)
    colloidal_aggregation, = aggregation_of_colloidal_iron(args.dissolved_organic_aggregation_parameter_1,
                                                           args.dissolved_organic_aggregation_parameter_2,
                                                           args.dissolved_organic_aggregation_parameter_3,
                                                           args.dissolved_organic_aggregation_parameter_4,
                                                           args.dissolved_organic_aggregation_parameter_5,
                                                           args.background_shear,
                                                           args.mixed_layer_shear,
                                                           args.z,
                                                           args.zₘₓₗ,
                                                           args.Fe,
                                                           Fe′,
                                                           args.DOC,
                                                           args.POC,
                                                           args.GOC)
    scavenging = iron_scavenging(λFe, args.POC + args.GOC, Fe′)
    bacterial_uptake = bacterial_iron_uptake(args.maximum_bacterial_growth_rate,
                                             args.particle_temperature_sensitivity,
                                             args.maximum_iron_ratio_in_bacteria,
                                             args.iron_half_saturation_for_bacteria,
                                             args.bacterial_iron_uptake_efficiency,
                                             args.T,
                                             args.Fe,
                                             Bact,
                                             LBact)

    ligand_aggregation_loss = ligand_aggregation(iron, args.Fe, args.DOC, args.T, λFe)

    return (small_particle_iron_remineralisation + grazing_waste + upper_trophic_waste -
            phytoplankton_iron_uptake - ligand_aggregation_loss - colloidal_aggregation -
            scavenging - bacterial_uptake)
end

@inline function free_iron(::SimpleIron, Fe, DOC, T)
    ligands = max(0.6, 0.09 * (DOC + 40) - 3)
    K = exp(16.27 - 1565.7 / max(T + 273.15, 5))
    Δ = 1 + K * ligands - K * Fe

    return (-Δ + √(Δ^2 + 4K * Fe)) / 2K
end

@inline function free_iron(iron::SimpleIron, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    DOC = @inbounds fields.DOC[i, j, k]
    Fe  = @inbounds fields.Fe[i, j, k]
    T   = @inbounds fields.T[i, j, k]

    return free_iron(iron, Fe, DOC, T)
end

end # module