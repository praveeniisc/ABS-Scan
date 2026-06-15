import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def setup_plot_style():
    """Configure modern, publication-quality plotting styles."""
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.titlesize": 16,
        "font.family": "sans-serif"
    })

def plot_results(csv_path, outdir):
    """
    Generate beautiful bar plots for both Binding Affinity Change (ddG) 
    and Structural Stability Change (dDOPE).
    """
    df = pd.read_csv(csv_path)
    
    # Filter out wild-type (WT) or receptor-only lines if present in the data
    mutants = df[df["Residue"] != "WT"].copy()
    if mutants.empty:
        return
        
    setup_plot_style()
    
    # 1. Plot Binding Affinity Contribution (ddG = G_mutant - G_WT)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create colors: Green for positive contribution (mutation weakens binding), Red for negative
    colors = ["#2ecc71" if val >= 0 else "#e74c3c" for val in mutants["ddG"]]
    
    sns.barplot(
        x="Residue",
        y="ddG",
        data=mutants,
        palette=colors,
        hue="Residue",
        legend=False,
        ax=ax,
        edgecolor="black",
        linewidth=0.8
    )
    
    ax.axhline(0, color="black", linestyle="-", linewidth=1.0)
    ax.set_ylabel(r"$\Delta\Delta G_{binding}$ (kcal/mol)", fontweight="bold")
    ax.set_xlabel("Binding Site Residues", fontweight="bold")
    ax.set_title(r"Residue Contribution to Binding Affinity ($\Delta\Delta G$)", pad=15, fontweight="bold")
    
    # Rotate residue labels if there are many of them
    if len(mutants) > 8:
        plt.xticks(rotation=45, ha="right")
        
    plt.tight_layout()
    affinity_plot_path = os.path.join(outdir, "binding_affinity_changes.png")
    plt.savefig(affinity_plot_path, dpi=300)
    plt.close()

    # 2. Plot DOPE Score Change (dDOPE = DOPE_mutant - DOPE_WT)
    if "dDOPE" in mutants.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Color: Red for destabilizing (positive dDOPE), Green for stabilizing (negative dDOPE)
        colors = ["#e74c3c" if val >= 0 else "#2ecc71" for val in mutants["dDOPE"]]
        
        sns.barplot(
            x="Residue",
            y="dDOPE",
            data=mutants,
            palette=colors,
            hue="Residue",
            legend=False,
            ax=ax,
            edgecolor="black",
            linewidth=0.8
        )
        
        ax.axhline(0, color="black", linestyle="-", linewidth=1.0)
        ax.set_ylabel(r"$\Delta$DOPE Score", fontweight="bold")
        ax.set_xlabel("Mutated Residues", fontweight="bold")
        ax.set_title(r"Effect of Alanine Mutation on Structural Stability ($\Delta$DOPE)", pad=15, fontweight="bold")
        
        if len(mutants) > 8:
            plt.xticks(rotation=45, ha="right")
            
        plt.tight_layout()
        dope_plot_path = os.path.join(outdir, "dope_score_changes.png")
        plt.savefig(dope_plot_path, dpi=300)
        plt.close()
