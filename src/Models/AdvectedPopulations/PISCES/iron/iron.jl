module Iron

export SimpleIron

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

@inline function iron_tendency(iron::SimpleIron,
                               minimum_iron_scavenging_rate,
                               load_specific_iron_scavenging_rate,
                               base_breakdown_rate,
                               particle_temperature_sensitivity,
                               maximum_iron_ratio_in_bacteria,
                               iron_half_saturation_for_bacteria,
                               bacterial_iron_uptake_efficiency,
                               maximum_bacterial_growth_rate,
                               dissolved_organic_aggregation_parameter_1,
                               dissolved_organic_aggregation_parameter_2,
                               dissolved_organic_aggregation_parameter_3,
                               dissolved_organic_aggregation_parameter_4,
                               dissolved_organic_aggregation_parameter_5,
                               microzooplankton_bacteria_concentration,
                               mesozooplankton_bacteria_concentration,
                               maximum_bacteria_concentration,
                               bacteria_concentration_depth_exponent,
                               doc_half_saturation_for_bacterial_activity,
                               nitrate_half_saturation_for_bacterial_activity,
                               ammonia_half_saturation_for_bacterial_activity,
                               phosphate_half_saturation_for_bacterial_activity,
                               iron_half_saturation_for_bacterial_activity,
                               nano_exudated_fraction,
                               nano_maximum_iron_ratio,
                               nano_half_saturation_for_iron_uptake,
                               nano_threshold_for_size_dependency,
                               nano_size_ratio,
                               nano_minimum_ammonium_half_saturation,
                               nano_minimum_nitrate_half_saturation,
                               nano_optimal_iron_quota,
                               nano_base_growth_rate,
                               nano_temperature_sensitivity,
                               diatom_exudated_fraction,
                               diatom_maximum_iron_ratio,
                               diatom_half_saturation_for_iron_uptake,
                               diatom_threshold_for_size_dependency,
                               diatom_size_ratio,
                               diatom_minimum_ammonium_half_saturation,
                               diatom_minimum_nitrate_half_saturation,
                               diatom_optimal_iron_quota,
                               diatom_base_growth_rate,
                               diatom_temperature_sensitivity,
                               zoo,
                               Fe,
                               DOC,
                               T,
                               POC,
                               GOC,
                               SFe,
                               CaCO₃,
                               PSi,
                               P,
                               PChl,
                               PFe,
                               D,
                               DChl,
                               DFe,
                               Z,
                               M,
                               NH₄,
                               NO₃,
                               PO₄,
                               Si,
                               O₂,
                               z,
                               zₘₓₗ,
                               zₑᵤ,
                               Si′,
                               background_shear,
                               mixed_layer_shear,
                               sinking_flux,
                               sinking_iron_flux,
                               first_anoxia_threshold,
                               second_anoxia_threshold)
    λFe = iron_scavenging_rate(minimum_iron_scavenging_rate,
                               load_specific_iron_scavenging_rate,
                               POC,
                               GOC,
                               CaCO₃,
                               PSi)
    Fe′ = free_iron(iron, Fe, DOC, T)

    food_availability = (; P, D, POC, Z)
    iron_availability = (; P = PFe / (P + eps(zero(P))),
                           D = DFe / (D + eps(zero(D))),
                           POC = SFe / (POC + eps(zero(POC))),
                           Z = zoo.micro.iron_ratio)

    Bact = bacteria_concentration(microzooplankton_bacteria_concentration,
                                 mesozooplankton_bacteria_concentration,
                                 maximum_bacteria_concentration,
                                 bacteria_concentration_depth_exponent,
                                 z,
                                 zₘₓₗ,
                                 zₑᵤ,
                                 Z,
                                 M)
    LBact = bacteria_activity(doc_half_saturation_for_bacterial_activity,
                              nitrate_half_saturation_for_bacterial_activity,
                              ammonia_half_saturation_for_bacterial_activity,
                              phosphate_half_saturation_for_bacterial_activity,
                              iron_half_saturation_for_bacterial_activity,
                              NH₄,
                              NO₃,
                              PO₄,
                              Fe,
                              DOC)

    O₂_min_1 = first_anoxia_threshold
    O₂_min_2 = second_anoxia_threshold
    
    small_particle_iron_remineralisation = degradation(Val(:SFe),
                                                        specific_degradation_rate(base_breakdown_rate,
                                                                                  particle_temperature_sensitivity,
                                                                                  O₂,
                                                                                  T,
                                                                                  O₂_min_1,
                                                                                  O₂_min_2),
                                                        SFe)
    grazing_waste = non_assimilated_iron(zoo, T, Z, M, food_availability, iron_availability, sinking_flux, sinking_iron_flux)
    upper_trophic_waste = upper_trophic_dissolved_iron(zoo, T, M)
    phytoplankton_iron_uptake = uptake(nano_exudated_fraction,
                                       nano_maximum_iron_ratio,
                                       nano_half_saturation_for_iron_uptake,
                                       nano_threshold_for_size_dependency,
                                       nano_size_ratio,
                                       nano_minimum_ammonium_half_saturation,
                                       nano_minimum_nitrate_half_saturation,
                                       nano_optimal_iron_quota,
                                       nano_base_growth_rate,
                                       nano_temperature_sensitivity,
                                       diatom_exudated_fraction,
                                       diatom_maximum_iron_ratio,
                                       diatom_half_saturation_for_iron_uptake,
                                       diatom_threshold_for_size_dependency,
                                       diatom_size_ratio,
                                       diatom_minimum_ammonium_half_saturation,
                                       diatom_minimum_nitrate_half_saturation,
                                       diatom_optimal_iron_quota,
                                       diatom_base_growth_rate,
                                       diatom_temperature_sensitivity,
                                       Val(:Fe),
                                       T,
                                       Fe,
                                       NO₃,
                                       NH₄,
                                       PO₄,
                                       Si,
                                       Si′,
                                       P,
                                       PChl,
                                       PFe,
                                       D,
                                       DChl,
                                       DFe)
    colloidal_aggregation, = aggregation_of_colloidal_iron(dissolved_organic_aggregation_parameter_1,
                                                           dissolved_organic_aggregation_parameter_2,
                                                           dissolved_organic_aggregation_parameter_3,
                                                           dissolved_organic_aggregation_parameter_4,
                                                           dissolved_organic_aggregation_parameter_5,
                                                           background_shear,
                                                           mixed_layer_shear,
                                                           z,
                                                           zₘₓₗ,
                                                           Fe,
                                                           Fe′,
                                                           DOC,
                                                           POC,
                                                           GOC)
    scavenging = iron_scavenging(λFe, POC + GOC, Fe′)
    bacterial_uptake = bacterial_iron_uptake(maximum_bacterial_growth_rate,
                                             particle_temperature_sensitivity,
                                             maximum_iron_ratio_in_bacteria,
                                             iron_half_saturation_for_bacteria,
                                             bacterial_iron_uptake_efficiency,
                                             T,
                                             Fe,
                                             Bact,
                                             LBact)

    ligand_aggregation_loss = ligand_aggregation(iron, Fe, DOC, T, λFe)

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