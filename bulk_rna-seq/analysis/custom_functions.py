import os
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np
import seaborn as sns
from adjustText import adjust_text
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from pydeseq2.ds import DeseqStats
from pydeseq2.dds import DeseqDataSet
import gseapy as gp
import textwrap

def dds_vst_function(data:pd.DataFrame=None, metadata:pd.DataFrame=None, condition:str=None):
    if condition is not None:
        data=data.loc[[c for c in data.index.values if c.split('_')[0] in condition],:]
        metadata=metadata.loc[[c for c in data.index.values if c.split('_')[0] in condition],:]
        
    #dds object
    dds= DeseqDataSet(counts=data.values,
                       metadata=metadata,
                       design_factors="Condition",
                       refit_cooks=True)
    dds.var_names=data.columns
    #vst function
    dds.vst()

    vst= pd.DataFrame(data=dds.layers['vst_counts'],
                      index=dds.obs.index,    
                      columns=data.columns.values)
    return dds, vst

def plot_transformations(data=None, vst=None,condition:str=None,save:bool=True, save_format:str='svg'):
    fig,ax = plt.subplots(1,3, figsize=(30,6))
    # non-normalized counts
    sns.scatterplot(x=data.mean(axis=0), y=data.var(axis=0), ax=ax[0])
    ax[0].set_ylim(-10e5,4*10e6)
    ax[0].set_xlim(-10e2, 40000)
    ax[0].set_ylabel('Gene variance')
    ax[0].set_xlabel('Raw Counts')
    ax[0].set_title('No Transformation')

    # vst transformed counts
    sns.regplot(x=vst.mean(axis=0), y=vst.var(axis=0),ax=ax[1],order=5, line_kws=dict(color="r"))
    ax[1].set_ylabel('Gene variance')
    ax[1].set_xlabel('Transformed counts')
    ax[1].set_title('VST')
    # ax[1].set_ylim(-0.1,0.5)

    # log1p transformed counts
    sns.regplot(x=np.log1p(data).mean(axis=0), y=np.log1p(data).var(axis=0),ax=ax[2],order=5, line_kws=dict(color="r"))
    ax[2].set_ylabel('Gene variance')
    ax[2].set_xlabel('Transformed counts')
    ax[2].set_title('Log1p Transformation')
    # ax[2].set_ylim(-0.1,3)

    if save:
        # save plot
        plt.savefig(f'results/figures/benchmarking_count_transformation_{condition}.{save_format}', dpi=300)
    plt.close()
    
    return fig


def distance_table(df1,df2,metric="euclidean"):
    df = pd.concat([df1,df2])
    x = squareform(pdist(df.values,metric=metric))
    x = np.nan_to_num(x)
    m = x.max()
    x = x[:df1.shape[0],:]
    x = x[:,df1.shape[0]:]

    return pd.DataFrame(x,index=df1.index.values,columns=df2.index.values), m

def plot_loadings_components(pca, component1, component2, features, n_genes, ax):
    """
    Function to plot the loadings of two specified PCA components on the same 2D plot
    and select the top genes based on the square mean of the loadings over all components.

    Parameters:
        pca (PCA): Fitted PCA object.
        component1 (int): Index of the first component to plot.
        component2 (int): Index of the second component to plot.
        features (list): List of feature names.
        n_genes (int): Number of top genes to select.

    Returns:
        None
    """
    # Calculate the mean square of loadings for each gene across all components
    loadings_sq_mean = np.mean(pca.components_ ** 2, axis=0)
    # Get the indices of the top genes based on the square mean
    top_gene_indices = np.argsort(loadings_sq_mean)[-n_genes:]

    # Get the loadings for the specified components
    loadings_comp1 = pca.components_[component1, top_gene_indices]
    loadings_comp2 = pca.components_[component2, top_gene_indices]

    # Calculate the minimum and maximum values of the loadings for both axes
    min_loading_x = min(loadings_comp1)
    max_loading_x = max(loadings_comp1)
    min_loading_y = min(loadings_comp2)
    max_loading_y = max(loadings_comp2)

    # Add a small buffer to the limits
    buffer_x = 0.1 * (max_loading_x - min_loading_x)
    min_loading_x -= buffer_x
    max_loading_x += buffer_x

    buffer_y = 0.1 * (max_loading_y - min_loading_y)
    min_loading_y -= buffer_y
    max_loading_y += buffer_y

    ax.set_xlim(min_loading_x, max_loading_x)  # Set x-axis limits
    ax.set_ylim(min_loading_y, max_loading_y)  # Set y-axis limits
    
    ax.quiver(np.zeros(n_genes), np.zeros(n_genes), loadings_comp1, loadings_comp2,
            angles='xy', scale_units='xy', scale=1, color='b', width=0.005, label=f'Components {component1+1}-{component2+1}')

    # Annotate the top genes based on loadings
    texts = []
    for idx in range(n_genes):
        x_pos = loadings_comp1[idx]
        y_pos = loadings_comp2[idx]
        text = features[top_gene_indices[idx]]

        texts.append(ax.text(x_pos, y_pos, text, fontsize=10, ha='left', va='bottom', color='black', bbox=dict(facecolor='white', alpha=0.8)))

    # Adjust text positions to avoid overlap with quiver vectors
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))

    ax.set_xlabel(f'Component {component1+1} Loading')
    ax.set_ylabel(f'Component {component2+1} Loading')
    ax.set_title(f'Top {n_genes} Genes based on Loadings of PCA Components {component1+1} and {component2+1}')
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.axvline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(False)
    ax.legend()

    return ax

def loadings_genelist(pca, component1:int=0, component2:int=1, features:np.ndarray=None, n_genes:int=100, path:str=None, format:str='csv'):
    """
    Function to save a table with the loadings of two specified PCA components based
    based on the square mean of the loadings over all components.

    Parameters:
        pca (PCA): Fitted PCA object.
        component1 (int): Index of the first component to plot.
        component2 (int): Index of the second component to plot.
        features (list): List of feature names.
        n_genes (int): Number of top genes to select.

    Returns:
        None
    """
    # empty df
    df = df_empty = pd.DataFrame({'Gene' : [], 'Loadings1': [], 'Loadings2': []})

    # Calculate the mean square of loadings for each gene across all components
    loadings_sq_mean = np.mean(pca.components_ ** 2, axis=0)

    # Get the indices of the top genes based on the square mean
    top_gene_indices = np.argsort(loadings_sq_mean)[-n_genes:]

    # Get the loadings for the specified components
    loadings_comp1 = pca.components_[component1, top_gene_indices]
    loadings_comp2 = pca.components_[component2, top_gene_indices]

    # construct df
    df.iloc[:,0] = features[top_gene_indices]
    df.iloc[:,1] = loadings_comp1
    df.iloc[:,2] = loadings_comp2

    df.to_csv(f'{path}/top{n_genes}_loadings.{format}')

    return

def plot_heatmap_markers(data:pd.DataFrame=None, samples=None, markers:list=None, palette:str='RdYlBu_r', 
                        figsize:tuple=(20,10), save:bool=True, filename:str=None, save_format:str='svg'):
    """
    Generates a heatmap with two subplots.

    Args:
        data (pd.DataFrame): The input data containing expression values.
        samples (list): List of sample names to filter the data.
        markers (list): List of marker genes to include in the heatmap.
        palette (str): Color palette for the heatmap.
        title (str): Title for the plot.
        figsize (tuple): width and height of the plot
        save (bool): Whether to save the figure.

    Returns:
        fig,ax
    """
    # check argument contents
    if data is not None:
        if markers is not None:
            markers=[m for m in markers if m in data.columns.values] 
            data=data.loc[:,markers]
        else:
            raise TypeError('Marker list not provided.')
        if samples is None:
            df=data 
        else:
            df=data.loc[[c for c in data.index.values if c[:-2] in samples],:]

        # mean for condition
        df_mean=df.groupby([x[:-2] for x in df.index.values]).mean().T
        df_mean=df_mean[samples]
    if data is None:
        raise TypeError('Count dataframe not provided.')

    # plot heatmap
    fig,ax=plt.subplots(1,2, figsize=(figsize[0],figsize[1]))

    # per sample heatmap
    sns.heatmap(df.T, ax=ax[0], cmap=palette)
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation = 90, fontsize = 12)
    ax[0].set_xlabel('')

    # mean expression heatmap
    sns.heatmap(df_mean, ax=ax[1],cmap=palette)
    ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation = 90, fontsize = 12)

    if save:
        plt.savefig(f'results/figures/heatmap_{filename}_markers.{save_format}', dpi=300)

def pairwise_comparison(dds, c1, c2, lfc_threshold: float = 1, padj_threshold: float = 0.05, save_format:str='csv'):
    """
    Performs a pairwise comparison between two conditions using DESeq2, and saves the results to CSV files.

    Parameters:
    - dds: DESeqDataSet object for DESeq2 analysis.
    - c1 (str): Name of the first condition to compare.
    - c2 (str): Name of the second condition to compare.
    - lfc_threshold (float): Threshold for log2 fold change to consider a gene as differentially expressed.
    - pvalue (float): P-value threshold to consider a gene as differentially expressed.

    Returns:
    - None
    """
    # Set up pair-wise comparison
    stat_res = DeseqStats(dds, contrast=('Condition', c1, c2), quiet=True)
    
    # Run statistical test: Wald test
    stat_res.summary()
    stat_res_df=stat_res.results_df
    # Get results dataframe and sort by log2 fold change (lfc)
    res = stat_res_df.sort_values(by='log2FoldChange', ascending=False)
    
    # Subset differentially expressed genes (DEGs)
    degs_df = res[(abs(res['log2FoldChange']) > lfc_threshold) & (res['padj'] < padj_threshold)]

    # Save full results and DEGs to CSV files
    res.to_csv(f"results/de/ranked_{c1}_vs_{c2}.{save_format}")
    degs_df.to_csv(f"results/de/degs_{c1}_vs_{c2}.{save_format}")

    return stat_res_df, degs_df



def plot_heatmap_degs(stat_res: pd.DataFrame = None, vst_counts: pd.DataFrame = None, c1: str = None, c2: str = None,
                      lfc_threshold: float = 1, pvalue_threshold: float = 0.05,
                      ntopgenes: int = 30, save_format: str = 'svg'):
    """
    Generates a heatmap for differentially expressed genes (DEGs) based on statistical results and variance-stabilized counts.

    Parameters:
    - stat_res (pd.DataFrame): DataFrame containing statistical results with columns 'log2FoldChange' and 'pvalue'.
    - vst_counts (pd.DataFrame): DataFrame of variance-stabilized transformed (VST) counts.
    - c1 (str): Name of the first condition to compare.
    - c2 (str): Name of the second condition to compare.
    - lfc_threshold (float): Threshold for log2 fold change to consider a gene as differentially expressed.
    - pvalue_threshold (float): P-value threshold to consider a gene as differentially expressed.
    - ntopgenes (int): Number of top genes to include in the heatmap.
    - save_format (str): File format for saving the heatmap (e.g., 'svg', 'png').

    Returns:
    - fig, ax: Matplotlib figure and axis objects of the generated heatmap.
    """
    # Sort by log2 fold change (lfc)
    res = stat_res.sort_values(by='log2FoldChange', ascending=False)
    
    # Subset differentially expressed genes (DEGs)
    degs_df = res[(abs(res['log2FoldChange']) > lfc_threshold) & (res['pvalue'] < pvalue_threshold)]

    ### Heatmap generation
    # Subset top DEGs by lfc, filtering out unwanted genes
    degs = degs_df.index.values
    degs = [d for d in degs if not d.startswith('Gm') and 'Rik' not in d]
    topgenes = np.concatenate([degs[:ntopgenes], degs[-ntopgenes:]])

    # Subset samples in comparisons and topgenes
    hdf = vst_counts.loc[[x for x in vst_counts.index.values if x[:-2] in [c1, c2]], topgenes]

    fig, ax = plt.subplots(1, 2, figsize=(20, 15))

    # Calculate mean expression for each condition
    hdf_mean = hdf.groupby(hdf.index.str[:-2]).mean().T
    hdf_mean = hdf_mean[[c1, c2]]
    
    # Plot heatmaps
    # Heatmap for individual replicates
    sns.heatmap(hdf.T, ax=ax[0], cmap='RdYlBu_r')
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=90, fontsize=12)
    ax[0].set_xlabel('')
    
    # Heatmap for mean expression
    sns.heatmap(hdf_mean, ax=ax[1], cmap='RdYlBu_r')
    ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=90, fontsize=12)
    ax[1].set_xlabel('')
    
    # Save heatmap figure
    plt.savefig(f'results/figures/heatmap_top{ntopgenes}_genes_{c1}_vs_{c2}.{save_format}')
    
    return fig, ax
    

def plot_volcano(res_df: pd.DataFrame = None, 
                 lfc_threshold: float = 1, padj_threshold: float = 0.05, basemean_filter: float = 10,
                markers: list = None, gene_labels:bool=True, size:float=20,
                palette1: list = ['red', 'green', 'blue', 'lightgrey'], palette2: list = ['purple'], title: str = '',
                get_degs:bool=False, save:bool=True, filename:str=None, save_format: str = 'svg', figsize:tuple=None):
    """
    Generates a volcano plot for visualizing differentially expressed genes (DEGs).

    Parameters:
    - res_df (pd.DataFrame): DataFrame containing the gene expression results with columns 'log2FoldChange', 'pvalue', and 'baseMean'.
    - lfc_threshold (float): Threshold for log2 fold change to consider a gene as differentially expressed.
    - pvalue_threshold (float): P-value threshold to consider a gene as differentially expressed.
    - basemean_filter (float): Minimum mean expression level to filter out lowly expressed genes.
    - markers (list): List of specific gene symbols to highlight on the plot.
    - gene_labels (bool): Whether to add gene labels to the plot.
    - size (float): Size for markers.
    - palette1 (list): List of colors for plotting categories: ['DEG', 'Log2FC', 'p-value', 'NS'].
    - palette2 (list): List of colors for highlighting specific markers.
    - title (str): Title of the plot.
    - get_degs (bool): Whether to save the list of differentially expressed genes (DEGs) to a CSV file.
    - save_format (str): File format for saving the plot (e.g., 'svg', 'png').

    Returns:
    - fig, ax: Matplotlib figure and axis objects of the generated plot.
    """
    
    # Check if res_df is provided
    if res_df is None:
        raise ValueError("A results DataFrame must be provided.")
        
    # Sorting the dataframe
    res = res_df.sort_values(by='log2FoldChange', ascending=False).copy()
    res['Symbol'] = res.index
    res['-log10padj'] = -np.log10(res['padj'])
    
    # Filter lowly expressed genes
    res = res[res['baseMean'] > basemean_filter]
    
    # Define color mapping function
    def map_color(row, lfc_threshold, padj_threshold, markers):
        if markers and row['Symbol'] in markers:
            return 'Marker'
        if abs(row['log2FoldChange']) > lfc_threshold and row['padj'] < padj_threshold:
            return 'DEG'
        if abs(row['log2FoldChange']) > lfc_threshold:
            return 'Log2FC'
        if row['padj'] < padj_threshold:
            return 'padj'
        return 'NS'
    
    # Apply color mapping
    res['Category'] = res.apply(map_color, axis=1, args=(lfc_threshold,padj_threshold, markers))
    
    # Plotting
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    if markers is not None:
        sns.scatterplot(data=res[res['Category'] == 'Marker'], x='log2FoldChange', y='-log10padj', hue='Category', s=size*2,
                        palette=palette2, edgecolor='black', ax=ax)

    
    sns.scatterplot(data=res, x='log2FoldChange', y='-log10padj', hue='Category', s=20,
                    hue_order=['DEG', 'Log2FC', 'padj', 'NS'],
                    palette=palette1, edgecolor='black', ax=ax)
    

    # Axes titles and lines
    ax.set_title(f'LFC threshold: {lfc_threshold} / padj threshold: {padj_threshold}', fontsize=8)
    ax.axhline(-np.log10(padj_threshold), zorder=0, color='black', linewidth=2, linestyle='--')
    ax.axvline(lfc_threshold, zorder=0, color='black', linewidth=2, linestyle='--')
    ax.axvline(-lfc_threshold, zorder=0, color='black', linewidth=2, linestyle='--')
    ax.legend(fontsize=8).set_title('')
    
    # Add gene labels
    if gene_labels:
        if markers is not None:
            if len(markers)>0:
                m_df = res[res['Category'] == 'Marker']
                labels = []
                for l in range(m_df.shape[0]):
                    x_pos = m_df['log2FoldChange'].iloc[l]
                    y_pos = m_df['-log10padj'].iloc[l]
                    label = m_df['Symbol'].iloc[l]
                    labels.append(ax.text(x_pos, y_pos, str(label), fontsize=8, ha='center', va='bottom', color='black', 
                                        bbox=dict(facecolor='white', alpha=0.5)))
                adjust_text(labels, arrowprops=dict(arrowstyle='-', color='gray'))
            else:
                raise ValueError('Empty marker list provided.')
    
    # Overall title
    plt.suptitle(title, y=0.96, x=0.52, fontsize=8)
    ax.set_xlabel('log2FoldChange',fontsize=8)
    ax.set_ylabel('-log10 p-adj value',fontsize=8)
    ax.tick_params(axis='y', which='major', labelsize=8)

    fig.tight_layout()
    # Save plot
    if save:
        plt.savefig(f'results/figures/volcano_{filename}.{save_format}', dpi=300,bbox_inches='tight')

    if get_degs:
        degs= res[(abs(res['log2FoldChange']) > lfc_threshold) & (res['padj'] < padj_threshold)]
        degs.to_csv(f'results/genelists/volcano_degs_{filename}.csv')

    return fig, ax

def go_enrich(f, lib, go_outdir):
    """
    Perform Gene Ontology (GO) enrichment analysis on differentially expressed genes (DEGs) and save results.

    Parameters:
    f (str): Path to the CSV file containing DEG data with log2 fold changes.
    lib (str): The gene set library to use for enrichment analysis.
    go_outdir (str): Directory where the enrichment results will be saved.

    Returns:
    enr_up: Enrichment results for upregulated genes.
    enr_down: Enrichment results for downregulated genes.
    """
    
    # Create the output directory for storing results if it does not exist
    os.makedirs(f'{go_outdir}/{lib}', exist_ok=True)
    
    # Read the DEG data from the CSV file, setting the first column as the index
    degs = pd.read_csv(f, index_col=0)

    # Extract the condition name from the file name (assuming the format 'degs_condition.csv')
    condition = f.split('degs_')[-1].split('.')[0]
    
    # Check if the DEG DataFrame is empty
    if degs.empty:
        print('No DEGs in ' + condition)
        return None, None  # Return None if there are no DEGs

    else:
        # Subset upregulated and downregulated genes based on log2 fold change
        degs_up = degs[degs.log2FoldChange > 0]
        degs_down = degs[degs.log2FoldChange < 0]

        # Perform GO enrichment analysis for upregulated genes
        enr_up = gp.enrichr(gene_list=degs_up.index.to_list(),
                            gene_sets=lib,
                            organism='mouse',
                            outdir=None, 
                            )
        # Save the enrichment results for upregulated genes to a CSV file
        enr_up.res2d.to_csv(f'{go_outdir}/{lib}/up_{condition}.csv', index=0)

        # Perform GO enrichment analysis for downregulated genes
        enr_down = gp.enrichr(gene_list=degs_down.index.to_list(),
                              gene_sets=lib,
                              organism='mouse',
                              outdir=None, 
                              )
        # Save the enrichment results for downregulated genes to a CSV file
        enr_down.res2d.to_csv(f'{go_outdir}/{lib}/down_{condition}.csv', index=0)

        # Return the enrichment results for both upregulated and downregulated genes
        return enr_up, enr_down


def plot_enrich(enr_obj, padj_threshold: float = 0.1, ntop:int=None, figsize: tuple = (15, 10),palette=None, title: str = None, filename: str = None, save_format: str = None):
    """
    Create and save enrichment plots based on the results of a Gene Set Enrichment Analysis.

    Parameters:
    enr_obj (object): An enrichment analysis object with results.
    padj_threshold (float): P-value threshold for significance (default is 0.1).
    figsize (tuple): Figure size (width, height) for the plots (default is (15, 10)).
    title (str): Title of the plot (default is None).
    filename (str): Base filename for saving the plots (default is None).
    save_format (str): Format for saving the plots (e.g., 'png', 'pdf') (default is None).

    Returns:
    None
    """
    
    # If no results are present, skip the plotting
    if enr_obj.results.shape[0] == 0:
        pass
    else:
        # Normalize p-values for color mapping
        norm = mpl.colors.Normalize(vmin=enr_obj.results['Adjusted P-value'].min(), vmax=enr_obj.results['Adjusted P-value'].max())
        cmapper = cm.ScalarMappable(norm=norm, cmap=palette)

        # Prepare data for barplot
        b_df = enr_obj.results.sort_values(by='Adjusted P-value', ascending=True)
        b_df.Term = b_df.Term.str.split(" \(GO").str[0]

        # Create a barplot with seaborn
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = sns.barplot(data=b_df.iloc[:ntop], x='Combined Score', y='Term', palette=cmapper.to_rgba(b_df['Adjusted P-value'].values), gap=0)
        ax.set_ylabel('')

        # Wrap long term names for readability
        ax.set_yticklabels([textwrap.fill(e, 40) for e in b_df['Term']], fontsize=8)
        plt.subplots_adjust(left=0.3)

        # Set the title and save the barplot
        plt.suptitle(title, y=0.90, fontsize=8)
        plt.savefig(f'{filename}_barplot.{save_format}', bbox_inches='tight')
        plt.close()