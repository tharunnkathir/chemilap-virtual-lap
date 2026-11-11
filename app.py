import streamlit as st

# Set page configuration (optional, but recommended)
st.set_page_config(page_title="ChemiLap Virtual Lab", layout="wide")

# Sidebar navigation menu
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Chemistry Quiz", "Virtual Lab", "About"])

# 1. The FIRST condition must ALWAYS start with 'if'
if page == "Chemistry Quiz":
    st.title("🧪 Chemistry Quiz")
    st.write("Welcome to the quiz section!")
    
    # Add your quiz code here
    answer = st.radio("What is the chemical symbol for Water?", ["H2O", "CO2", "O2"])
    if st.button("Submit Answer"):
        if answer == "H2O":
            st.success("Correct! 🎉")
        else:
            st.error("Incorrect, try again.")

# 2. Subsequent conditions use 'elif'
elif page == "Virtual Lab":
    st.title("🔬 Virtual Lab")
    st.write("Perform your chemical experiments here.")
    
    # Add your virtual lab code here

elif page == "About":
    st.title("ℹ️ About ChemiLap")
    st.write("This app is designed to help students learn chemistry interactively.")
