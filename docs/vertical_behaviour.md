# Vertical behaviour of *Pelagia noctiluca*: evidence and model design

Where a jellyfish is in the water column decides which current carries it. In the Alboran Sea the
surface Atlantic Jet can flow at 0.5–1 m s⁻¹ while the water 150 m below moves slowly or even in
another direction, so an animal that spends the day at depth and the night at the surface does not
follow the surface drift. This note summarises what is known about the vertical behaviour of
*P. noctiluca*, and how (and how not) it is represented in `jellydrift`.

## 1. Evidence

| Topic | Evidence | Source |
|---|---|---|
| Diel vertical migration (DVM) | Adults are near the surface at night and deep by day. In the Ligurian Sea they were rarely seen at the surface in daylight but always present offshore at night. | Ferraris et al. (2012); Berline et al. (2013) |
| Day depth | Below 300 m in the NW Mediterranean; 185 m (New England) and 90 m (Florida) from submersibles. | Franqueville (1971) and Larson et al. (1991), as cited in Berline et al. (2013) |
| Day depth (recent, in situ) | ROV footage from 660 Mediterranean sites: live animals mostly between 108 and 190 m; all individuals observed below ~195 m were dead. | Toma et al. (2026) |
| Night depth | Particles held at 5 m at night in a Ligurian Sea stranding model. | Berline et al. (2013) |
| Vertical swimming speed | ~2 m min⁻¹ (≈120 m h⁻¹); a 5→305 m migration modelled in 3 h (100 m h⁻¹). An IBM uses 1800–2400 m day⁻¹ (75–100 m h⁻¹) of active swimming for adults, plus a small negative buoyancy. | Arai (1997) as cited in Berline et al. (2013); El Rahi et al. (2020, model code) |
| Timing | Clock-based in Berline et al. (at depth 09:00–21:00, at 5 m 23:00–06:00); synchronised with sunrise/sunset in the IBMs. | Berline et al. (2013); El Rahi et al. (2020); Vodopivec & Malej (2025) |
| Shallow water | Over bottoms shallower than 305 m animals were kept 100 m above the seabed; over bottoms shallower than 100 m they stop migrating and stay at 5 m, because bottom contact and wave stress disturb migration. Another IBM caps the dive depth as a function of bottom depth (≈25 m over 100 m, 100 m over 300 m). | Berline et al. (2013); El Rahi et al. (2020, model code) |
| Early life stages | Ephyrae stay in the upper 50 m (99.9 %): mostly in the upper 25 m at night, spread over 0–50 m by day, never below the thermocline. Transport models of early stages treat them as surface particles. | Pastor-Prieto et al. (2021); Bergamasco et al. (2022) |
| Temperature | Pulsation rate increases linearly with temperature between 8 and 22 °C and peaks at 16–19 °C; outside 8–22 °C active swarming is not possible (Gulf of Trieste). An IBM scales swimming speed with temperature. | Avian et al. (1991); El Rahi et al. (2020, model code) |
| Waves, wind, turbulence | On the Moroccan side of the Alboran Sea, strandings correlate with wave height and wind direction (n = 33). No study found that documents *P. noctiluca* diving to avoid rough seas. | Aouititen et al. (2019) |
| Moonlight | No study found for *P. noctiluca*. | — |

## 2. What the model does

`jellydrift.behaviour.DielMigration` (parameters in `config/default.yaml`, section `behaviour`):

1. **Light cue.** Each individual heads down when the sun is above the horizon at its position and
   up when it is below (sun-synchronised, as in the IBMs, rather than fixed clock hours).
2. **Depths.** Night target 5 m. The day target is a preferred depth drawn per individual from a
   normal distribution (mean = `day_depth_m`, relative s.d. 0.2, truncated to [5 m, 2 × mean]),
   representing the spread of observed depths. The mean is a **sweep parameter**: 50 m (a shallow,
   ephyra-like migration), 150 m (the depth band of live animals in the ROV survey) and 300 m (the
   classic value used by Berline et al.), plus a passive surface-drifter control.
3. **Finite speed.** Individuals move towards the target at most 100 m h⁻¹, so dawn descents and
   dusk ascents take one to three hours and cross the intermediate current layers.
4. **Seabed.** Individuals stay at least 100 m above the seabed (model bathymetry) and do not
   migrate at all where the bottom is shallower than 100 m (Berline et al., 2013).
5. **Currents at depth.** Horizontal currents are interpolated linearly in depth between the
   hourly surface field and daily fields at the model levels nearest to 10, 25, 50, 100, 150, 200
   and 300 m (10.5, 26.2, 51.4, 97.9, 153.4, 203.2 and 303.6 m).
6. **Stranding only from the upper layer.** An individual can only strand when it is in the upper
   25 m (the night layer). Deeper individuals whose step would reach land, which happens where
   the coast is steep and a 4 km coastal cell is already deep, keep their position for that step.

## 3. What the model does not do, and why

- **Temperature-dependent swimming.** In the reanalysis for 2022 (basin means over
  5.3°W–2°W, 35.5°N–36.6°N), water at 100–300 m stays between 13.6 and 15.8 °C all year, and
  monthly mean surface temperature ranges from 15.3 °C (February) to 24.5 °C (August), with
  daily extremes of 14.2–28.7 °C. Day depths are therefore inside the 8–22 °C active range,
  but from July to September the night layer is warmer than 22 °C, where Avian et al. (1991)
  report that active swimming is impaired. Since temperature mainly acts on the vertical
  speed, that speed is tested directly in a sensitivity experiment (50 and 200 m h⁻¹) instead
  of adding 3-D temperature fields. Summer impairment remains an open question for future work.
- **Life stages.** Only adults, the stinging stage people see on beaches, are represented.
  Ephyrae would behave roughly like the 50 m scenario.
- **Storm avoidance and moonlight.** No species-specific evidence was found, so neither is
  modelled.
- **Wind and waves.** These are horizontal drivers rather than behaviour, and remain a
  limitation: wave-induced Stokes drift and direct wind drag on surface individuals are not
  included (see the report's Limitations).

## References

- Aouititen, M., Bekkali, R., Nachit, D., Luan, X., & Mrhraoui, M. (2019). Predicting jellyfish
  strandings in the Moroccan north-west Mediterranean coastline. *European Scientific Journal*,
  15(2), 72. https://doi.org/10.19044/esj.2019.v15n2p72
- Arai, M. N. (1997). *A Functional Biology of Scyphozoa*. Chapman & Hall, London.
- Avian, M., Rottini Sandrini, L., & Stravisi, F. (1991). The effect of seawater temperature on the
  swimming activity of *Pelagia noctiluca* (Forsskål). *Bollettino di Zoologia*, 58(2), 135–143.
  https://doi.org/10.1080/11250009109355744
- Bergamasco, A., Cucco, A., Guglielmo, L., et al. (2022). Observing and modeling long-term
  persistence of *P. noctiluca* in coupled complementary marine systems (Southern Tyrrhenian Sea
  and Messina Strait). *Scientific Reports*, 12, 14905. https://doi.org/10.1038/s41598-022-18832-2
- Berline, L., Zakardjian, B., Molcard, A., Ourmières, Y., & Guihou, K. (2013). Modeling jellyfish
  *Pelagia noctiluca* transport and stranding in the Ligurian Sea. *Marine Pollution Bulletin*,
  70(1–2), 90–99. https://doi.org/10.1016/j.marpolbul.2013.02.016
- El Rahi, J., Weeber, M. P., & El Serafy, G. (2020). Modelling the effect of behavior on the
  distribution of the jellyfish Mauve stinger (*Pelagia noctiluca*) in the Balearic Sea using an
  individual-based model. *Ecological Modelling*, 433, 109230.
  https://doi.org/10.1016/j.ecolmodel.2020.109230 — model code:
  https://github.com/joeelrahi/Pelagianoctiluca (parameters read from the code; no code reused).
- Ferraris, M., Berline, L., Lombard, F., Guidi, L., Elineau, A., Mendoza-Vera, J. M., Lilley,
  M. K. S., Taillandier, V., & Gorsky, G. (2012). Distribution of *Pelagia noctiluca* (Cnidaria,
  Scyphozoa) in the Ligurian Sea (NW Mediterranean Sea). *Journal of Plankton Research*, 34(10),
  874–885. https://doi.org/10.1093/plankt/fbs049
- Franqueville, C. (1971). Macroplancton profond (invertébrés) de la Méditerranée nord-occidentale.
  *Tethys*, 3, 11–56.
- Larson, R. J., Mills, C. E., & Harbison, G. R. (1991). Western Atlantic midwater hydrozoan and
  scyphozoan medusae: in situ studies using manned submersibles. *Hydrobiologia*, 216/217, 311–317.
- Pastor-Prieto, M., Bahamon, N., Sabatés, A., Canepa, A., Gili, J.-M., Carreton, M., & Company,
  J. B. (2021). Spatial heterogeneity of *Pelagia noctiluca* ephyrae linked to water masses in the
  Western Mediterranean. *PLoS ONE*, 16(4), e0249756. https://doi.org/10.1371/journal.pone.0249756
- Toma, M., Rivela, A., Angiolillo, M., Bavestrello, G., Canese, S., Cau, A., & Bo, M. (2026).
  Insights into the presence and ecology of *Pelagia noctiluca* (Forsskål, 1775) in the deep
  Mediterranean Sea. *Marine Environmental Research*, 221, 108328.
  https://doi.org/10.1016/j.marenvres.2026.108328
- Vodopivec, M., & Malej, A. (2025). Spatially explicit individual-based model reveals the mauve
  stinger jellyfish distribution in the Mediterranean Sea. *Ecological Modelling*, 505, 111109.
  https://doi.org/10.1016/j.ecolmodel.2025.111109
