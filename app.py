elif page == "Chemistry Quiz":

    st.title("📝 Chemistry Knowledge Quiz")

    questions = [
        {
            "question": "What is the chemical formula of water?",
            "options": ["CO2", "H2O", "O2", "CH4"],
            "answer": "H2O"
        },
        {
            "question": "What is the molecular formula of benzene?",
            "options": ["C6H6", "CH4", "C2H6O", "H2O"],
            "answer": "C6H6"
        },
        {
            "question": "Which compound is commonly found in vinegar?",
            "options": ["Methane", "Benzene", "Acetic Acid", "Glucose"],
            "answer": "Acetic Acid"
        },
        {
            "question": "What does SMILES represent?",
            "options": [
                "A programming language",
                "A molecular representation system",
                "A chemical element",
                "A laboratory instrument"
            ],
            "answer": "A molecular representation system"
        },
        {
            "question": "What is the chemical formula of methane?",
            "options": ["CH4", "CO2", "C6H6", "H2O"],
            "answer": "CH4"
        },
        {
            "question": "Which element is represented by the symbol O?",
            "options": ["Gold", "Oxygen", "Silver", "Iron"],
            "answer": "Oxygen"
        },
        {
            "question": "What is the molecular formula of glucose?",
            "options": ["C6H12O6", "H2O", "CH4", "C6H6"],
            "answer": "C6H12O6"
        },
        {
            "question": "Which software library is used for cheminformatics in this project?",
            "options": ["NumPy", "RDKit", "Matplotlib", "Flask"],
            "answer": "RDKit"
        },
        {
            "question": "What does molecular weight measure?",
            "options": [
                "Molecular mass",
                "Temperature",
                "Pressure",
                "Volume"
            ],
            "answer": "Molecular mass"
        },
        {
            "question": "What is the formula of ethanol?",
            "options": ["C2H6O", "C6H6", "CH4", "CO2"],
            "answer": "C2H6O"
        }
    ]

    with st.form("quiz_form"):

        answers = []

        for i, q in enumerate(questions):

            st.write(f"**Q{i+1}. {q['question']}**")

            answer = st.radio(
                "Choose your answer:",
                q["options"],
                key=f"question_{i}"
            )

            answers.append(answer)

        submitted = st.form_submit_button("Submit Quiz")

    if submitted:

        score = sum(
            answers[i] == questions[i]["answer"]
            for i in range(len(questions))
        )

        st.success(f"Your Score: {score}/10")

        if score >= 8:
            st.balloons()
            st.write("Excellent work!")

        elif score >= 5:
            st.write("Good effort! Keep learning.")

        else:
            st.write("Keep practising and try again.")

        with st.expander("View Correct Answers"):

            for i, q in enumerate(questions):
                st.write(
                    f"Q{i+1}: Correct Answer: {q['answer']}"
                )
