"""
    SimpleIron(; excess_scavenging_enhancement = 1000)

Parameterisation for iron evolution, not the "complex chemistry" model
of Aumount et al, 2015. Iron is scavenged (i.e. permanently removed from
the model) when the free iron concentration exceeds the ligand concentration
at a rate modified by `excess_scavenging_enhancement`.
"""
@kwdef struct SimpleIron{FT}
    excess_scavenging_enhancement :: FT = 1000.0 # unitless
     maximum_ligand_concentration :: FT = 0.6    # μmol Fe / m³
           dissolved_ligand_ratio :: FT = 0.09   # μmol Fe / mmol C
end

required_biogeochemical_tracers(::SimpleIron) = tuple(:Fe)

const SimpleIronPISCES = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:SimpleIron}

@inline function iron_inputs(bgc::SimpleIronPISCES, i, j, k, grid, clock, fields, auxiliary_fields)
    Fe = @inbounds fields.Fe[i, j, k]
    DOC = @inbounds fields.DOC[i, j, k]

    T = @inbounds fields.T[i, j, k]

    scavenging_rate = iron_scavenging_rate(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    colloidal_aggregation, = aggregation_of_colloidal_iron(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    scavenging = iron_scavenging(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    bacterial_uptake = bacterial_iron_uptake(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    small_particle_iron_remineralisation = degredation(bgc.particulate_organic_matter, Val(:SFe), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    phytoplankton_iron_uptake = uptake(bgc.phytoplankton, Val(:Fe), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    grazing_waste = non_assimilated_iron(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    upper_trophic_waste = upper_trophic_dissolved_iron(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    Lₜᵐᵃˣ = bgc.iron.maximum_ligand_concentration
    Lₜ = bgc.iron.dissolved_ligand_ratio * DOC - Lₜᵐᵃˣ
    total_ligand_concentration = max(Lₜᵐᵃˣ, Lₜ)

    free_iron_concentration = free_iron(bgc.iron, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    excess_iron = max(zero(Fe), Fe - total_ligand_concentration)

    ligand_aggregation_loss = bgc.iron.excess_scavenging_enhancement *
                              scavenging_rate *
                              excess_iron *
                              free_iron_concentration

    return IronInputs(
        small_particle_iron_remineralisation,
        grazing_waste,
        upper_trophic_waste,
        phytoplankton_iron_uptake,
        ligand_aggregation_loss,
        colloidal_aggregation,
        scavenging,
        bacterial_uptake,
    )
end

@inline function (bgc::SimpleIronPISCES)(i, j, k, grid, ::Val{:Fe}, clock, fields, auxiliary_fields)
    inputs = iron_inputs(bgc, i, j, k, grid, clock, fields, auxiliary_fields)

    return iron_tendency(bgc.iron, inputs)
end