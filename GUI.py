# ----------------------------------------------------- #
# Welcome the code for the GUI of the project. This     #
# code opens a GUI where the user chooses the reactions #
# he desires and what parameters he wishes to use.      # 
# ----------------------------------------------------- # 


# Import of Libraries
import tkinter as tk
from tkinter import messagebox
import subprocess
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#Choice of colors for plots
colors_colourblind = np.array(["blue", "black", "orange", "cyan", "palevioletred", "lime", "darkmagenta"])


# ------------------------------------------------ #
# Functions for Compilation, Running, and Plotting #
# ------------------------------------------------ #


def compile_cpp_code(cpp_file):

    """
    This function serves to compile the C++ code to run the simulation of
    Plasma-Surface Recombination. It assumes the user has the header files
    in a folder named inc2 and the .cpp files are in a folder named src2 .

    This function takes a .cpp file and outputs a executable.
    """

    try:

        compile_command = f"g++ -I inc2 src2/*.cpp -o exec"
        result = subprocess.run(
            compile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        if result.returncode != 0:
            messagebox.showerror("Compilation Error", f"Error compiling code:\n{result.stderr.decode()}")
            return None
        
        else:
            return "./exec"
        
    except Exception as e:

        messagebox.showerror("Error", f"Error compiling code: {str(e)}")
        return None

def run_cpp_code(react_program, selected_reactions, parameters):

    """
    This function serves to run the C++ code compiled in the previous function.

    This function takes a executable, a string of the selected reactions
    and of the parameters and outputs None. 
    """

    try:

        command = [react_program] + selected_reactions + parameters

        process = subprocess.Popen(command, stderr=subprocess.PIPE, text=True, bufsize=1)

        process.wait()
        error_text = process.stderr.read()

        if process.returncode != 0:
            messagebox.showerror("Error", f"Error running code:\n{error_text}")

        else:
            messagebox.showinfo("Success", "Code executed successfully. Check .txt files for the output")
            output_file = "output.txt"
            if os.path.exists(output_file):
                generate_plots(output_file)
            else:
                print("Output file not found; plot generation skipped.")

    except Exception as e:

        messagebox.showerror("Error", f"Error: {str(e)}")



def generate_plots(file_path):

    """
    This function serves to generate plots regarding the chosen reactions
    by the user. The functions are made to be, by default, normalized. The
    concentrations evolution is plotted here.

    This function takes a .cpp file and outputs a executable.
    """

    try:

        data = pd.read_csv(file_path, delimiter='\t')
        time = data.iloc[:, 0]

        pop_cols = [col for col in data.columns if col.startswith("Population")]
        rate_cols = [col for col in data.columns if col.startswith("R")]

        
        plt.figure(figsize=(9, 16))
        for i, col in enumerate(pop_cols):
            plt.plot(time, data[col]/np.max(data[col]), label=col, color=colors_colourblind[i % len(colors_colourblind)])
        
        plt.xlabel('Time [s]', fontsize=18)
        plt.ylabel('Concentration', fontsize=18)
        plt.xscale('log')
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlim(10e-16, 10e-11)

        plt.legend(fontsize=14)
        plt.grid(True)
        plt.show()


        plt.figure(figsize=(9, 16))
        for i, col in enumerate(rate_cols):
            plt.plot(time, data[col]/np.max(data[col]), label=col, color=colors_colourblind[i % len(colors_colourblind)])
        
        plt.xlabel('Time [s]', fontsize=18)
        plt.ylabel('Reaction Rate', fontsize=18)
        plt.xscale('log')
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlim(10e-16, 10e-11)

        plt.legend(fontsize=14)
        plt.grid(True)
        plt.show()

    except Exception as e:

        messagebox.showerror("Error", f"Error generating plots: {str(e)}")

# These will be the reactions chosen
selected_reactions_global = []


# -------------------------------------------------- #
# Dictionaries that maps each type of reaction and   #
# their concentration, their reaction parameters and # 
# the stop time.                                     #
# -------------------------------------------------- #


reaction_parameters = {

    "Basic": [
        {"label": "Constant of A → B", "default": "0.1"},
        {"label": "Constant of B → A", "default": "0.05"},
        {"label": "A Concentration", "default": "10000"},
        {"label": "B Concentration", "default": "5000"},
        {"label": "Stop Time", "default": "80.0"}
    ],

    "Physisorption": [
        {"label": "General Parameter Tw", "default": "200"},
        {"label": "General Parameter Tg", "default": "500"},
        {"label": "General Parameter M", "default": "16e-3"},
        {"label": "Sticking Probability for Physisorption", "default": "1.0"},
        {"label": "Desorption Frequency", "default": "10e15"},
        {"label": "Desorption Barrier", "default": "30e3"},
        {"label": "A Concentration", "default": "1e5"},
        {"label": "Fᵥ Concentration", "default": "1.5e5"},
        {"label": "Aₚ Concentration", "default": "0.0"},
        {"label": "Stop Time", "default": "5e-13"}
    ],

    "Chemisorption": [
        {"label": "General Parameter Tw", "default": "200"},
        {"label": "General Parameter Tg", "default": "500"},
        {"label": "General Parameter M", "default": "16e-3"},
        {"label": "Sticking Probability for Chemisorption", "default": "1.0"},
        {"label": "Pre-Exp. Factor for Recomb.", "default": "1.0"},
        {"label": "Recomb. Barrier", "default": "17.5e3"},
        {"label": "A Concentration", "default": "1e5"},
        {"label": "Sᵥ Concentration", "default": "3e3"},
        {"label": "Aₛ Concentration", "default": "0.0"},
        {"label": "A₂ Concentration", "default": "0.0"},
        {"label": "Stop Time", "default": "5e-13"}
    ],

    "Surface Diffusion": [
        {"label": "General Parameter Tw", "default": "200"},
        {"label": "General Parameter Tg", "default": "500"},
        {"label": "General Parameter M", "default": "16e-3"},
        {"label": "Diffusion Frequency", "default": "10e13"},
        {"label": "Diffusion Barrier", "default": "15e3"},
        {"label": "Aₚ Concentration", "default": "0.0"},
        {"label": "Aₛ Concentration", "default": "0.0"},
        {"label": "Fᵥ Concentration", "default": "1.5e5"},
        {"label": "Sᵥ Concentration", "default": "3e3"},
        {"label": "Stop Time", "default": "5e-13"}
    ],

    "Langmuir-Hinshelwood recombination": [
        {"label": "General Parameter Tw", "default": "200"},
        {"label": "General Parameter Tg", "default": "500"},
        {"label": "General Parameter M", "default": "16e-3"},
        {"label": "Diffusion Frequency", "default": "10e13"},
        {"label": "Diffusion Barrier", "default": "15e3"},
        {"label": "Pre-Exp. Factor for Recomb.", "default": "1.0"},
        {"label": "Recomb. Barrier", "default": "17.5e3"},
        {"label": "Recomb. Barrier for Two Physisorbed Atoms", "default": "0.0"},
        {"label": "Aₚ Concentration", "default": "0.0"},
        {"label": "Aₛ Concentration", "default": "0.0"},
        {"label": "Fᵥ Concentration", "default": "1.5e5"},
        {"label": "Sᵥ Concentration", "default": "3e3"},
        {"label": "A₂ Concentration", "default": "0.0"},
        {"label": "Stop Time", "default": "5e-13"}
    ]
}

# Fixed mapping for concentration labels to species names.
label_to_species = {
    "A Concentration": "A",
    "B Concentration": "B",
    "Aₚ Concentration": "Af",
    "Fᵥ Concentration": "Fv",
    "Aₛ Concentration": "As",
    "Sᵥ Concentration": "Sv",
    "A₂ Concentration": "A2"
}

# The reactions might have the same type of species. As
# such, one must define the order of each type.
fixedOrder = ["A", "B", "Af", "As", "Fv", "Sv", "A2"]


# ------------------------------------ #
# First GUI Window: Reaction Selection #
# ------------------------------------ #


def show_reaction_selection_window():

    """
    This function shows the first window where the user will
    select the reactions he wishes to use.
    """

    reaction_win = tk.Tk()
    reaction_win.title("Select Reactions")
    tk.Label(reaction_win, text="Select Reactions:").pack(pady=10)
    
    reaction_vars = {}
    reaction_choices = list(reaction_parameters.keys())
    for choice in reaction_choices:
        var = tk.BooleanVar()
        reaction_vars[choice] = var
        tk.Checkbutton(reaction_win, text=choice, variable=var).pack(pady=2, anchor="w")
    
    def next_step():
        global selected_reactions_global
        selected = [r for r, var in reaction_vars.items() if var.get()]

        # The Basic reaction is just an example and, as such, it does not
        # make sense to be used in conjunction with other reactions.
        if "Basic" in selected and len(selected) > 1:
            messagebox.showerror("Error", "Basic reaction cannot be mixed with any other reaction.")
            return
        
        if not selected:
            messagebox.showerror("Error", "Please select at least one reaction")

        else:
            selected_reactions_global = selected
            reaction_win.destroy()
            show_parameter_window()
    
    tk.Button(reaction_win, text="Next", command=next_step).pack(pady=10)
    reaction_win.mainloop()


# --------------------------------------- #
# Second GUI Window: Choice of Parameters #
# --------------------------------------- #


def show_parameter_window():

    """
    This function shows the second window where the user will
    select the parameters he wishes to use for the chosen reactions.
    """

    param_win = tk.Tk()
    param_win.title("Enter Simulation Parameters")
    param_win.geometry("400x500")

    # Create a canvas.
    canvas = tk.Canvas(param_win)
    scrollbar = tk.Scrollbar(param_win, orient="vertical", command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)
    scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")
    
    # Grouping parameters
    general_params = {}
    constants = {}
    stop_time_param = None
    
    # The window has different types of parameters seperated
    for reaction in selected_reactions_global:
        for param in reaction_parameters[reaction]:
            label = param["label"]

            if "Stop Time" in label:
                if stop_time_param is None:
                    stop_time_param = param

            elif "Concentration" in label:
                continue  

            elif "General Parameter" in label:
                if label not in general_params:
                    general_params[label] = param

            else:
                if label not in constants:
                    constants[label] = param

    union_species = set()
    for reaction in selected_reactions_global:
        for param in reaction_parameters[reaction]:
            if "Concentration" in param["label"]:
                if param["label"] in label_to_species:
                    union_species.add(label_to_species[param["label"]])

    union_species_ordered = [s for s in fixedOrder if s in union_species]
    
    species_defaults = {}
    for reaction in selected_reactions_global:
        for param in reaction_parameters[reaction]:
            if "Concentration" in param["label"]:
                species = label_to_species.get(param["label"])
                if species and species not in species_defaults:
                    species_defaults[species] = param["default"]

    entry_widgets = []

    def create_group_frame(group_name, param_list):
        if not param_list:
            return
        frame = tk.LabelFrame(scrollable_frame, text=group_name, padx=10, pady=10)
        frame.pack(padx=10, pady=5, fill="x")
        for param in param_list:
            tk.Label(frame, text=param["label"] + ":").pack(pady=2, anchor="w")
            entry = tk.Entry(frame)
            entry.insert(0, param["default"])
            entry.pack(pady=2, anchor="w", fill="x")
            entry_widgets.append(entry)

    # Parameters that are general and do not depend on the reactions chosen
    create_group_frame("General Parameters", list(general_params.values()))
    
    # Parameters specific to the reactions.
    create_group_frame("Reaction Constants", list(constants.values()))
    
    # Concentration of each species.
    if union_species_ordered:
        frame = tk.LabelFrame(scrollable_frame, text="Concentrations", padx=10, pady=10)
        frame.pack(padx=10, pady=5, fill="x")
        for species in union_species_ordered:
            tk.Label(frame, text=f"{species} Concentration:").pack(pady=2, anchor="w")
            entry = tk.Entry(frame)
            default_val = species_defaults.get(species, "")
            entry.insert(0, default_val)
            entry.pack(pady=2, anchor="w", fill="x")
            entry_widgets.append(entry)
    
    # Stop time of the reaction.
    if stop_time_param is not None:
        frame = tk.LabelFrame(scrollable_frame, text="Stop Time", padx=10, pady=10)
        frame.pack(padx=10, pady=5, fill="x")
        tk.Label(frame, text=stop_time_param["label"] + ":").pack(pady=2, anchor="w")
        entry = tk.Entry(frame)
        entry.insert(0, stop_time_param["default"])
        entry.pack(pady=2, anchor="w", fill="x")
        entry_widgets.append(entry)
    
    # Function to call the function that runs the simulation
    def on_run():
        cpp_file = "src2/Plasma-Surface-Recombination.cpp"
        compiled_program = compile_cpp_code(cpp_file)

        if not compiled_program:
            messagebox.showerror("Error", "Compilation failed.")
            return

        all_values = [e.get().strip() for e in entry_widgets]
        run_cpp_code(compiled_program, selected_reactions_global, all_values)
    
    tk.Button(scrollable_frame, text="Run Code", command=on_run).pack(pady=10)
    param_win.mainloop()

# Start the GUI
show_reaction_selection_window()
