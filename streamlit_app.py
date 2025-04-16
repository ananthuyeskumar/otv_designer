from otv_flight_dynamics import *

st.set_page_config(layout="wide")
st.title("OTV Designer")

# Sidebar
st.sidebar.title("App Modules")
module = st.sidebar.radio("Select a module", ["Flight Dynamics", "OTV Sizing", "Commercial"])

if module == "Flight Dynamics":
    st.header("Flight Dynamics Requirements")
    delta_v_budget = 0

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Altitude Change")
        h1 = st.number_input("Enter Initial Altitude (km)", 300, 40000, value=500)
        h2 = st.number_input("Enter Final Altitude (km)", 300, 40000, value=600)
        result = choose_transfer(h1, h2)
        st.markdown(f"<span style='color:red; font-size:16px;'>**Delta-V Required (km/s)**: {result['delta_v_kms']:.5f}</span>", unsafe_allow_html=True)
        delta_v_budget += result['delta_v_kms']

        st.subheader("Inclination Change")
        h_incl = st.number_input("Enter Altitude (km)", 300, 40000, value=500)
        inc = st.number_input("Enter Change in Inclination (deg)", 0.0, 30.0, value=5.0)
        dv_inc = inclination_change_delta_v(h_incl, inc)
        st.markdown(f"<span style='color:red; font-size:16px;'>**Delta-V Required (km/s)**: {dv_inc:.5f}</span>", unsafe_allow_html=True)
        delta_v_budget += dv_inc

    with col2:
        st.subheader("LTAN Change")
        h_ltan = st.number_input("Enter Initial Altitude (km)", 300, 40000, value=500, key='h_ltan')
        ltan_i = st.number_input("Enter Initial LTAN (hours)", 0.0, 24.0, value=10.5, key='ltan_i')
        ltan_f = st.number_input("Enter Final LTAN (hours)", 0.0, 24.0, value=6.0, key='ltan_f')
        ltan_dur = st.number_input("Enter duration (days)", 10, 360, value=90, key='ltan_dur')

        result = ltan_drift_adjustment(
            initial_ltan_hours=ltan_i,
            target_ltan_hours=ltan_f,
            initial_altitude_km=h_ltan,
            duration_days=ltan_dur
        )
        st.markdown(f"<span style='color:red; font-size:16px;'>**Delta-V Required (km/s)**: {result['total_delta_v_kms']:.5f}</span>", unsafe_allow_html=True)
        delta_v_budget += result['total_delta_v_kms']
        st.markdown(f"<span style='color:blue; font-size:16px;'>**Drift Orbit Inclination (deg)**: {result['new_inclination_deg']:.2f}</span>", unsafe_allow_html=True)
        st.markdown(f"<span style='color:blue; font-size:16px;'>**Drift Orbit Altitude (km)**: {result['new_altitude_km']:.2f}</span>", unsafe_allow_html=True)

    
    st.subheader("Phase Maneuver")
    ph_col1, ph_col2 = st.columns(2)
    with ph_col1:
        h_phase = st.number_input("Initial Orbit Altitude (km)", 300, 2000, value=450, key="h_phase")
        phase_angle = st.number_input("Desired Phase Angle (deg)", 1.0, 359.0, value=90.0)
    with ph_col2:
        num_orbits = st.number_input("Number of Orbits to Phase", 1, 100, value=18)

    phase_result = find_phasing_orbit(
        h_initial_km=h_phase,
        phase_deg=phase_angle,
        num_orbits=num_orbits
    )

    st.markdown(f"<span style='color:red; font-size:16px;'>**Delta-V Required (km/s)**: {phase_result['total_delta_v_km_s']:.5f}</span>", unsafe_allow_html=True)
    st.markdown(f"<span style='color:blue; font-size:16px;'>**Phasing Orbit Altitude (km)**: {phase_result['phasing_orbit_altitude_km']:.2f}</span>", unsafe_allow_html=True)
    st.markdown(f"<span style='color:blue; font-size:16px;'>**Drift Orbit Period (min)**: {phase_result['phasing_orbit_period_min']:.2f}</span>", unsafe_allow_html=True)
    st.markdown(f"<span style='color:blue; font-size:16px;'>**Total Drift Time (min)**: {phase_result['total_drift_time_min']:.2f}</span>", unsafe_allow_html=True)
    st.markdown(f"<span style='color:purple; font-size:16px;'>**Phasing Orbit Apoapsis (km)**: {phase_result['apoapsis_km']:.2f}</span>", unsafe_allow_html=True)
    st.markdown(f"<span style='color:purple; font-size:16px;'>**Phasing Orbit Periapsis (km)**: {phase_result['periapsis_km']:.2f}</span>", unsafe_allow_html=True)
    st.markdown(f"<span style='color:orange; font-size:15px;'>Note: The phasing orbit is elliptical, chosen to match the desired phase angle over {num_orbits} orbits.</span>", unsafe_allow_html=True)

    delta_v_budget += phase_result['total_delta_v_km_s']


    # Summary
    st.markdown("---")
    st.metric(label="Total Delta-V Required (km/s):", value=f"{delta_v_budget:.5f}")
