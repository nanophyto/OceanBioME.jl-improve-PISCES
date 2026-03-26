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

@inline function iron_tendency(bgc::SimpleIronPISCES, i, j, k, grid, clock, fields, auxiliary_fields)
    Fe = @inbounds fields.Fe[i, j, k]
    DOC = @inbounds fields.DOC[i, j, k]
    T = @inbounds fields.T[i, j, k]
    O₂ = @inbounds fields.O₂[i, j, k]

    POC = @inbounds fields.POC[i, j, k]
    GOC = @inbounds fields.GOC[i, j, k]
    SFe = @inbounds fields.SFe[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]
    PSi = @inbounds fields.PSi[i, j, k]

    P = @inbounds fields.P[i, j, k]
    PChl = @inbounds fields.PChl[i, j, k]
    PFe = @inbounds fields.PFe[i, j, k]
    D = @inbounds fields.D[i, j, k]
    DChl = @inbounds fields.DChl[i, j, k]
    DFe = @inbounds fields.DFe[i, j, k]

    Z = @inbounds fields.Z[i, j, k]
    M = @inbounds fields.M[i, j, k]
    NH₄ = @inbounds fields.NH₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    PO₄ = @inbounds fields.PO₄[i, j, k]
    Si = @inbounds fields.Si[i, j, k]

    zₘₓₗ = @inbounds auxiliary_fields.zₘₓₗ[i, j, k]
    zₑᵤ = @inbounds auxiliary_fields.zₑᵤ[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())
    Si′ = @inbounds bgc.silicate_climatology[i, j, k]

    sinking_flux = edible_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)
    sinking_iron_flux = edible_iron_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    pom = bgc.particulate_organic_matter

    return iron_tendency(bgc.iron,
                         pom.minimum_iron_scavenging_rate,
                         pom.load_specific_iron_scavenging_rate,
                         pom.base_breakdown_rate,
                         pom.temperature_sensitivity,
                         pom.maximum_iron_ratio_in_bacteria,
                         pom.iron_half_saturation_for_bacteria,
                         pom.bacterial_iron_uptake_efficiency,
                         pom.maximum_bacterial_growth_rate,
                         bgc.dissolved_organic_matter,
                         bgc.phytoplankton,
                         bgc.zooplankton,
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
                         bgc.background_shear,
                         bgc.mixed_layer_shear,
                         sinking_flux,
                         sinking_iron_flux,
                         bgc.first_anoxia_threshold,
                         bgc.second_anoxia_threshold)
end

@inline function (bgc::SimpleIronPISCES)(i, j, k, grid, ::Val{:Fe}, clock, fields, auxiliary_fields)
    return iron_tendency(bgc, i, j, k, grid, clock, fields, auxiliary_fields)
end
