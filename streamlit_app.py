from otv_flight_dynamics import *

st.title("OTV Configuration Designer")


st.header("Fligth Dynamics Requirements")

st.subheader("Altitude Change")

h1 = st.number_input("Enter Initial Altitude", value=500)
h2 = st.number_input("Enter Final Altitude", value=600)

delta_v_budget = 0

result = choose_transfer(h1, h2)

st.write("Delta-V Required", result['delta_v_kms'])

delta_v_budget += result['delta_v_kms']

st.subheader("Inclination Change")

h_incl = st.number_input("Enter Altitude", value=500)
inc = st.number_input("Enter Change in Inclination", value=5)

dv_inc = inclination_change_delta_v(h_incl, inc)  

st.write("Delta-V Required", dv_inc)

delta_v_budget += dv_inc

st.subheader("LTAN Change")

h_ltan = st.number_input("Enter Initial Altitude", value=500, key='h_ltan')
ltan_i = st.number_input("Enter Initial LTAN (hours)", value=9, key='ltan_i')
ltan_f = st.number_input("Enter Final LTAN (hours)", value=6, key='ltan_f')
ltan_dur = st.number_input("Enter duration (days)", value=90, key= 'ltan_dur')

result = ltan_drift_adjustment(
    initial_ltan_hours=ltan_i,
    target_ltan_hours=ltan_f,
    initial_altitude_km=h_ltan,
    duration_days=ltan_dur
)


st.write("Delta-V Required", result['total_delta_v_kms'])

delta_v_budget += result['total_delta_v_kms']

st.write("Total Delta-V Required", delta_v_budget)