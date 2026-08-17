import streamlit as st

st.set_page_config(page_title="Dwell-time analysis tools", layout="centered")

pages = [
    st.Page("views/home.py",          title="Home", default=True),
    st.Page("views/cm_model.py",      title="CM model"),
    st.Page("views/dm_model.py",      title="DM model"),
    st.Page("views/cm_simulation.py", title="CM simulation"),
    st.Page("views/dm_simulation.py", title="DM simulation"),
    st.Page("views/cm_3state.py",     title="CM 3-state"),
    st.Page("views/dm_3state.py",     title="DM 3-state"),
]
st.navigation(pages).run()
