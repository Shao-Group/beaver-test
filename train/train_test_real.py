# %%
import csv
import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import make_scorer, precision_score, recall_score, precision_recall_curve, accuracy_score, roc_auc_score
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, cross_validate, StratifiedKFold
import ast
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import dump
import statistics
from sklearn.utils import resample
import joblib
import PyPDF2
import math
import os

# Constants
BEAVER = 'beaver-test'
HUMAN_DIR = '../results/HEK293T_results/'
HUMAN_BEAVER_DIR = HUMAN_DIR + BEAVER
HUMAN_SAMPLE_SIZE_LIST = [2, 5, 10, 20, 30, 50, 80, 100, 150, 192]
#HUMAN_SAMPLE_SIZE_LIST = [2]
MOUSE_DIR = '../results/Mouse-Fibroblast_results/'
MOUSE_BEAVER_DIR = MOUSE_DIR + BEAVER
MOUSE_SAMPLE_SIZE_LIST = [2, 5, 10, 30, 50, 80, 100, 150, 200, 300, 369]
#MOUSE_SAMPLE_SIZE_LIST = [2]

THRESHOLD = 0.5
TRAIN_BASE_DIR = f'../train'
TRAIN_DIR = f'{TRAIN_BASE_DIR}/realDataset-{BEAVER}'
TRAIN_CHRM = ['1','2','3','4','5','6','7','8','9']
os.makedirs(TRAIN_DIR, exist_ok=True)

HUMAN_SAMPLE_SIZE_PLOT_LIST = [5, 10, 30, 50, 100, 192]
MOUSE_SAMPLE_SIZE_PLOT_LIST = [10, 30, 50, 100, 200, 369]

# %%
def process_tmap_files_other_chrs_noSpecificGT(tool, sample_size, base_dir, species):
    """
    Process .tmap files for a specific tool and sample size, and generate a performance summary.
    """
    tool_paths = {
        'transmeta': lambda cell_id: f"{base_dir}/transmeta/merge_{sample_size}/individual.{cell_id}.other.TransMeta.bam{cell_id}.otherchrm.gtf.tmap",
        'psiclass': lambda cell_id: f"{base_dir}/psiclass/merge_{sample_size}/cov_1/individual.{cell_id}.other.psiclass_sample_{cell_id}.otherchrm.gtf.tmap",
        'aletsch': lambda cell_id: f"{base_dir}/aletsch/merge_{sample_size}/gtf/{cell_id}.other.{cell_id}.otherchrm.gtf.tmap",
        'stringtie2': lambda cell_id: f"{base_dir}/stringtie2/{cell_id}.other.{cell_id}.otherchrm.gtf.tmap",
        'scallop2': lambda cell_id: f"{base_dir}/scallop2/{cell_id}.other.{cell_id}.otherchrm.gtf.tmap"
    }
    
    path_function = tool_paths.get(tool)
    if not path_function:
        raise ValueError(f"Unsupported tool: {tool}")
    
    # Determine the range of cell_ids based on the tool
    start_id = 0 if tool in ['psiclass', 'aletsch'] else 1
    end_id = sample_size - 1 if tool in ['psiclass', 'aletsch'] else sample_size
    
    all_data = []
    cell_stats = []
    for cell_id in range(start_id, end_id + 1):
        tmap_file = path_function(cell_id)
        gt_file_id = cell_id + 1 if tool in ['psiclass', 'aletsch'] else cell_id
        #gt_file = f"{HUMAN_SPECIFIC_GT}/{gt_file_id}.list" if species == 'HEK293T' else f"{MOUSE_SPECIFIC_GT}/{gt_file_id}.list"
        
        if not os.path.exists(tmap_file):
            print(f"Warning: File not found - {tmap_file}")
            continue
        
        # Load cell-specific ground truth
        #with open(gt_file, 'r') as f:
            #cell_specific_gt = set(line.strip() for line in f)
        
        # Process tmap file
        cell_data = []
        with open(tmap_file, 'r') as tmap:
            next(tmap)  # skip the header
            for line in tmap:
                fields = line.strip().split('\t')
                ref_gene_id, ref_id, class_code, qry_id = fields[0], fields[1], fields[2], fields[4]
                
                specific_gt = True #ref_id in cell_specific_gt
                label_general = class_code == '='
                label_specific = label_general and label_general
                
                cell_data.append({
                    'data_source': species,
                    'sample_size': sample_size,
                    'cell_id': gt_file_id,
                    'ref_gene_id': ref_gene_id,
                    'ref_id': ref_id,
                    'class_code': class_code,
                    'qry_id': qry_id,
                    'specificGT': specific_gt,
                    'label_general': label_general,
                    'label_specific': label_specific
                })
        
        # Calculate statistics for this cell
        total = len(cell_data)
        matching_general = sum(d['label_general'] for d in cell_data)
        matching_specific = sum(d['label_specific'] for d in cell_data)
        precision_general = matching_general / total if total > 0 else 0
        precision_specific = matching_specific / total if total > 0 else 0
        
        cell_stats.append({
            'data_source': species,
            'sample_size': sample_size,
            'cell_id': gt_file_id,
            'total': total,
            'matching_general': matching_general,
            'precision_general': precision_general,
            'matching_specific': matching_specific,
            'precision_specific': precision_specific
        })

        print(cell_stats[-1])
        
        all_data.extend(cell_data)
        
        # Save detailed tmap file output
        output_dir = f"{TRAIN_DIR}/tools-other-chrs/{tool}/{species}"
        os.makedirs(output_dir, exist_ok=True)
        output_file = f"{output_dir}/{species}.{sample_size}.cell-{gt_file_id}.csv"
        pd.DataFrame(cell_data).to_csv(output_file, index=False)
        #print(f"Saved detailed tmap output for {tool}, {species}, sample size {sample_size}, cell {gt_file_id} to {output_file}")
    
    return all_data, cell_stats

def process_all_tools_other_chrs(species):
    """
    Process all tools and sample sizes, calculate statistics, and save results to CSV files.
    """
    base_dir = HUMAN_DIR if species == 'HEK293T' else MOUSE_DIR
    sample_sizes = HUMAN_SAMPLE_SIZE_LIST if species == 'HEK293T' else MOUSE_SAMPLE_SIZE_LIST
    tools = ['transmeta', 'psiclass', 'aletsch', 'stringtie2', 'scallop2']
    #tools = ['aletsch-link']
    
    for tool in tools:
        all_stats = []
        for sample_size in sample_sizes:
            _, cell_stats = process_tmap_files_other_chrs_noSpecificGT(tool, sample_size, base_dir, species)
            all_stats.extend(cell_stats)
        
        # Save statistics to CSV
        output_dir = f"{TRAIN_BASE_DIR}/realDataset-tools-otherchrs/{tool}"
        os.makedirs(output_dir, exist_ok=True)
        output_file = f"{output_dir}/{species}_statistics.csv"
        pd.DataFrame(all_stats).to_csv(output_file, index=False)
        print(f"Saved statistics for {tool}, {species} to {output_file}")

# %%
def calc_stats(values):
    """Calculate statistics for a list of values."""
    if not values:
        return [0, 0, 0, 0, 0]  # Default values if the list is empty
    return [round(x, 3) for x in [np.min(values), np.max(values), np.median(values), np.mean(values), np.std(values)]]


def calc_stats_nonZero(values):
    """Calculate statistics for non-zero values in a list."""
    non_zero_values = [v for v in values if v != 0]
    
    if not non_zero_values:
        return [0, 0, 0, 0, 0]  # Default values if there are no non-zero values
    
    return [
        round(np.min(non_zero_values), 3),
        round(np.max(non_zero_values), 3),
        round(np.median(non_zero_values), 3),
        round(np.mean(non_zero_values), 3),
        round(np.std(non_zero_values), 3) if len(non_zero_values) > 1 else 0
    ]

def get_meta_feature_vector(v, sample_size):
    """Extract feature vector from raw data."""
    features = []
    
    # Basic information
    transcript_id, chrom = v[0], v[1]
    features.extend([transcript_id, chrom])
    
    # Transcript support information
    bottleneck_coverage = round(float(v[2]), 3)
    highest_coverage = round(float(v[3]), 3)
    extendable_score = round(float(v[4]), 3)
    num_junctions = int(v[5])
    features.extend([bottleneck_coverage, highest_coverage, extendable_score, num_junctions])
    
    # Compatible junction coverage
    comp_junc_cov = [round(float(x) / math.log(float(sample_size)), 3) for x in v[6:6+num_junctions]]
    features.extend(calc_stats(comp_junc_cov))
    
    # Junction coverage (regardless of compatibility)
    junc_cov = [round(float(x) / math.log(float(sample_size)), 3) for x in v[6+num_junctions:6+2*num_junctions]]
    features.extend(calc_stats(junc_cov))
    
    # Cell support information
    cells_support = int(v[6+2*num_junctions])
    if(cells_support == 0): return None
    features.append(cells_support)
    features.append(round(cells_support / float(sample_size), 3))
    
    # Cell support details
    cell_support_details = v[7+2*num_junctions:7+2*num_junctions+2*cells_support]
    cell_junction_counts = [int(cell_support_details[i+1]) for i in range(0, len(cell_support_details), 2)]
    features.extend(calc_stats(cell_junction_counts))

    # Number of cells covering all junctions
    cells_covering_all_junctions = sum(1 for count in cell_junction_counts if count == num_junctions)
    features.append(round(cells_covering_all_junctions / float(sample_size), 3))
    features.append(round(cells_covering_all_junctions / float(cells_support), 3))

    # Fragment information
    num_fragments = int(v[7+2*num_junctions+2*cells_support])
    features.append(num_fragments)
    
    # Fragment coverage
    fragment_coverage = [round(float(x), 3) for x in v[8+2*num_junctions+2*cells_support:]]
    features.extend(calc_stats(fragment_coverage))
    
    return features

def get_labels_from_tmap(tmap_file):
    """Extract labels from tmap file."""
    labels = {}
    with open(tmap_file, 'r') as tmap:
        next(tmap)  # skip the header
        for line in tmap:
            fields = line.strip().split('\t')
            qry_id, class_code = fields[4], fields[2]
            labels[qry_id] = 1 if class_code == '=' else 0
    return labels

def process_meta_data(data_dir, data_source, sample_size, train_chrm):
    """Process data for a specific sample size."""
    all_data = []
    
    tmap_file = f'{data_dir}/merge_{sample_size}/merge_{sample_size}.linkmerge_{sample_size}.gtf.tmap'
    feature_file = f'{data_dir}/merge_{sample_size}/linkmerge_{sample_size}_feature.csv'
    labels = get_labels_from_tmap(tmap_file)
        
    # Process feature file
    with open(feature_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            features = get_meta_feature_vector(row, sample_size)
            if(features == None): continue
            transcript_id, chrom = features[0], features[1]
            if transcript_id in labels:
                data_row = [data_source, sample_size] + features + [labels[transcript_id]]
                all_data.append(data_row)
            else:
                print(f"Warning: Inconsistency between feature and tmap at {transcript_id}")
    
    # Create DataFrame
    feature_names = [
        'data_source', 'sample_size', 'transcript_id', 'chr',
        'bottleneck_coverage', 'highest_coverage', 'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        #'min_cell_junc_norm', 'max_cell_junc_norm', 'median_cell_junc_norm', 'mean_cell_junc_norm', 'std_cell_junc_norm',
        #'most_occ_junc_num', 'most_occ_junc_num_freq',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov',
        'label_general'
    ]
    df = pd.DataFrame(all_data, columns=feature_names)
    
    # Split data into train and test sets
    train_df = df[df['chr'].isin(train_chrm)]
    test_df = df[~df['chr'].isin(train_chrm)]
    
    return train_df, test_df

def process_and_save_meta_data(data_dir, data_source, sample_size_list, train_chrm):
    """Process and save data for all sample sizes."""
    for sample_size in sample_size_list:
        train_df, test_df = process_meta_data(data_dir, data_source, sample_size, train_chrm)
        
        os.makedirs(f'{TRAIN_DIR}/train_test_meta_data', exist_ok=True)
        
        train_df.to_csv(f'{TRAIN_DIR}/train_test_meta_data/{data_source}_{sample_size}_train.csv', index=False)
        test_df.to_csv(f'{TRAIN_DIR}/train_test_meta_data/{data_source}_{sample_size}_test.csv', index=False)
        
        print(f"Saved {data_source} data for sample size {sample_size}")

# %%
def load_and_combine_meta_data(train_dir, data_sources):
    """Load, print statistics, and combine train and test data from CSV files."""
    train_data = []
    test_data = []
    
    # Define column types, specifying 'chr' as string
    dtype_dict = {'chr': str}
    
    print("\nData Statistics:")
    
    for source in data_sources:
        if source == 'HEK293T':
            sample_sizes = HUMAN_SAMPLE_SIZE_LIST
        elif source == 'Mouse-Fibroblast':
            sample_sizes = MOUSE_SAMPLE_SIZE_LIST
        else:
            raise ValueError(f"Unknown data source: {source}")
        
        for size in sample_sizes:
            train_file = f'{train_dir}/train_test_meta_data/{source}_{size}_train.csv'
            test_file = f'{train_dir}/train_test_meta_data/{source}_{size}_test.csv'
            
            train_df = pd.read_csv(train_file, dtype=dtype_dict)
            test_df = pd.read_csv(test_file, dtype=dtype_dict)
            
            train_data.append(train_df)
            test_data.append(test_df)
            
            print(f"\n{source} - Sample size {size}:")
            print(f"Train - Positive: {train_df['label_general'].sum()}, Negative: {len(train_df) - train_df['label_general'].sum()}")
            print(f"Test - Positive: {test_df['label_general'].sum()}, Negative: {len(test_df) - test_df['label_general'].sum()}")
            print(f"Initial Precision: {train_df['label_general'].sum()+test_df['label_general'].sum()}/{len(train_df)+len(test_df)}")
    combined_train = pd.concat(train_data, ignore_index=True)
    combined_test = pd.concat(test_data, ignore_index=True)
    
    print(f"\nTotal train data shape: {combined_train.shape}, positive: {combined_train['label_general'].sum()}, negative: {combined_train.shape[0] - combined_train['label_general'].sum()}")
    print(f"Total test data shape: {combined_test.shape}, positive: {combined_test['label_general'].sum()}, negative: {combined_test.shape[0] - combined_test['label_general'].sum()}")
    
    return combined_train, combined_test
    
def train_and_evaluate_meta_model(X_train, y_train, X_test, y_test, model_name):
    """Train a Random Forest model, evaluate its performance, and attach probabilities to test data."""
    feature_columns = [
        'sample_size',
        'bottleneck_coverage', 'highest_coverage', 'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        #'min_cell_junc_norm', 'max_cell_junc_norm', 'median_cell_junc_norm', 'mean_cell_junc_norm', 'std_cell_junc_norm',
        #'most_occ_junc_num', 'most_occ_junc_num_freq',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov'
    ]
    
    X_train = X_train[feature_columns]
    X_test = X_test[feature_columns]

    class_weights = {0: 1, 1: 5}
    model = RandomForestClassifier(n_estimators=100, max_depth=12, random_state=42, class_weight=class_weights)
    model.fit(X_train, y_train)
    
    # Train performance
    y_train_pred = model.predict(X_train)
    y_train_prob = model.predict_proba(X_train)[:, 1]
    train_precision = precision_score(y_train, y_train_pred)
    train_recall = recall_score(y_train, y_train_pred)
    train_roc_auc = roc_auc_score(y_train, y_train_prob)
    
    print(f'{model_name} Train Results:')
    print(f'Train Precision: {train_precision:.3f}')
    print(f'Train Recall: {train_recall:.3f}')
    print(f'Train ROC curve (area = {train_roc_auc:.3f})')
    
    # Test performance
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    test_precision = precision_score(y_test, y_pred)
    test_recall = recall_score(y_test, y_pred)
    test_roc_auc = roc_auc_score(y_test, y_prob)
    
    print(f'\n{model_name} Test Results:')
    print(f'Test Precision: {test_precision:.3f}')
    print(f'Test Recall: {test_recall:.3f}')
    print(f'Test ROC curve (area = {test_roc_auc:.3f})')
    
    # Save the model
    model_filename = f'../models/{model_name}_roc={test_roc_auc:.3f}.joblib'
    joblib.dump(model, model_filename)
    
    return model, model_filename, y_prob, y_train_prob

# %%
def process_meta_individual_cells(combined_test, data_dir, data_source, sample_size):
    """
    Process individual cell data using precomputed predictions and save results.
    """
    cell_dir = f"{data_dir}/merge_{sample_size}/linkmerge_{sample_size}_sgtf"
    output_dir = f"{TRAIN_DIR}/predictions-meta/{data_source}/{sample_size}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter combined_test for the current data_source and sample_size
    current_test = combined_test[(combined_test['data_source'] == data_source) & 
                                     (combined_test['sample_size'] == sample_size)]
    
    for cell_id in range(1, sample_size + 1):
        tmap_file = os.path.join(cell_dir, f"individual.{cell_id}.{cell_id}.gtf.tmap")
        output_file = os.path.join(output_dir, f"{data_source}.{sample_size}.cell-{cell_id}.csv")
        
        if not os.path.exists(tmap_file):
            print(f"Warning: tmap file for cell {cell_id} in {data_source}, sample size {sample_size} does not exist.")
            continue
        
        # Read tmap file
        tmap_data = []
        with open(tmap_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                fields = line.strip().split('\t')
                ref_gene_id, ref_id, class_code, qry_id = fields[0], fields[1], fields[2], fields[4]
                label = 1 if class_code == '=' else 0
                tmap_data.append([ref_gene_id, ref_id, class_code, qry_id, label])
        
        tmap_df = pd.DataFrame(tmap_data, columns=['ref_gene_id', 'ref_id', 'class_code', 'qry_id', 'label'])
        
        # Get intersection of transcript_ids
        cell_transcripts = set(tmap_df['qry_id'])
        test_transcripts = set(current_test['transcript_id'])
        common_transcripts = cell_transcripts.intersection(test_transcripts)
        
        # Filter tmap_df and current_test for common transcripts
        tmap_df = tmap_df[tmap_df['qry_id'].isin(common_transcripts)]
        cell_test = current_test[current_test['transcript_id'].isin(common_transcripts)]
        
        # Merge predictions with tmap data
        merged_df = pd.merge(tmap_df, cell_test[['transcript_id', 'y_prob_general', 'label_general']], 
                             left_on='qry_id', right_on='transcript_id')
        
        # Check label consistency
        label_mismatch = (merged_df['label'] != merged_df['label_general']).sum()
        if label_mismatch > 0:
            print(f"Warning: {label_mismatch} label mismatches found in cell {cell_id}.")
            merged_df = merged_df[merged_df['label'] == merged_df['label_general']]
        
        # Sort by y_prob in descending order
        merged_df = merged_df.sort_values('y_prob_general', ascending=False).reset_index(drop=True)
        
        # Calculate cumulative matches and precision
        merged_df['match_cum_general'] = merged_df['label_general'].cumsum()
        merged_df['precision_cum_general'] = merged_df['match_cum_general'] / (merged_df.index + 1)
        
        # Select and order columns
        result_df = merged_df[['ref_gene_id', 'ref_id', 'class_code', 'qry_id', 'label_general', 'y_prob_general', 'match_cum_general', 'precision_cum_general']]
        
        # Save results
        result_df.to_csv(output_file, index=False)
        
        # Print statistics
        total_instances = len(result_df)
        positive_cases = result_df['label_general'].sum()
        negative_cases = total_instances - positive_cases
        print(f"Cell {cell_id} in {data_source}, sample size {sample_size}:")
        print(f"  Total instances: {total_instances}")
        #print(f"  Positive cases: {positive_cases}")
        #print(f"  Negative cases: {negative_cases}")

# %%
def count_non_zero(values):
    """Count non-zero values in a list."""
    return sum(1 for v in values if v != 0)

def continuous_non_zero(values):
    """Calculate max and min continuous non-zero values."""
    if not values:
        return 0, 0
    
    counts = []
    current_count = 0
    
    for v in values:
        if v != 0:
            current_count += 1
        else:
            if current_count > 0:
                counts.append(current_count)
                current_count = 0
    
    if current_count > 0:
        counts.append(current_count)
    
    if not counts:
        return 0, 0
    
    return max(counts), min(counts)
    
def get_specific_feature_vector(v, sample_size):
    """Extract feature vector from raw data."""
    features = []
    
    # Basic information
    transcript_id, chrom = v[0], v[1]
    features.extend([transcript_id, chrom])
    
    # Input and output counts
    input_genes, input_transcripts, output_transcripts = int(v[2]), int(v[3]), int(v[4])
    features.extend([input_genes, input_transcripts, output_transcripts])
    
    # Transcript support information
    supporting_junctions = int(v[5])
    ratio_supporting_junctions = round(float(v[6]), 3)
    features.extend([supporting_junctions, ratio_supporting_junctions])
    # Number of junctions
    num_junctions = int(v[7])
    #features.append(num_junctions)
    
    # Process cell_comp_junc_cov
    start_idx = 8
    end_idx = start_idx + num_junctions
    cell_comp_junc_cov = [float(x) for x in v[start_idx:end_idx]]
    features.extend(calc_stats(cell_comp_junc_cov))
    
    # New features for cell_comp_junc_cov
    features.append(count_non_zero(cell_comp_junc_cov))
    max_continuous, min_continuous = continuous_non_zero(cell_comp_junc_cov)
    features.extend([max_continuous, min_continuous])
    
    # Process cell_junc_cov
    start_idx = end_idx
    end_idx = start_idx + num_junctions
    cell_junc_cov = [float(x) for x in v[start_idx:end_idx]]
    features.extend(calc_stats(cell_junc_cov))
    
    # New features for cell_junc_cov
    features.append(count_non_zero(cell_junc_cov))
    max_continuous, min_continuous = continuous_non_zero(cell_junc_cov)
    features.extend([max_continuous, min_continuous])
    
    return features

# %%
def get_specific_labels_from_tmap_noSpecificGT(tmap_file):
    """Extract labels from tmap file."""
    labels = {}
    with open(tmap_file, 'r') as tmap:
        next(tmap)  # skip the header
        for line in tmap:
            fields = line.strip().split('\t')
            ref_id, class_code, qry_id = fields[1], fields[2], fields[4]
            label_general = class_code == '='
            label_specific = label_general
            labels[qry_id] = {
                'label_general': label_general,
                'label_specific': label_specific
            }
    return labels

def combine_general_to_specific_noSpecificGT(data_dir, data_source, sample_size, train_chrm, meta_prediction):
    """Process data for a specific sample size."""
    all_data = []
    
    for cell_id in range(1, sample_size + 1):
        tmap_file = f"{data_dir}/merge_{sample_size}/linkmerge_{sample_size}_sgtf/individual.{cell_id}.{cell_id}.gtf.tmap"
        feature_file = f"{data_dir}/merge_{sample_size}/linkmerge_{sample_size}_sgtf/{cell_id}_feature.csv"
        #gt_file = f"{HUMAN_SPECIFIC_GT}/{cell_id}.list" if data_source == 'HEK293T' else f"{MOUSE_SPECIFIC_GT}/{cell_id}.list"
        
        if not os.path.exists(tmap_file) or not os.path.exists(feature_file): #or not os.path.exists(gt_file):
            print(f"Warning: File not found for cell {cell_id}")
            continue
        
        # Load cell-specific ground truth
        #with open(gt_file, 'r') as f:
            #cell_specific_gt = set(line.strip() for line in f)
        
        # Get labels from tmap file
        labels = get_specific_labels_from_tmap_noSpecificGT(tmap_file)
        
        # Process feature file
        with open(feature_file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                features = get_specific_feature_vector(row, sample_size)
                if features is None:
                    continue
                
                transcript_id, chrom = features[0], features[1]
                if transcript_id in labels:
                    data_row = [data_source, sample_size, cell_id, transcript_id, chrom] + features[2:] + [
                        labels[transcript_id]['label_general'],
                        labels[transcript_id]['label_specific']
                    ]
                    all_data.append(data_row)
    
    # Create DataFrame
    feature_names = [
        'data_source', 'sample_size', 'cell_id', 'transcript_id', 'chr',
        'input_genes', 'input_transcripts', 'output_transcripts',
        'cell_supporting_junctions', 'cell_ratio_supporting_junctions',
        'cell_min_comp_cov', 'cell_max_comp_cov', 'cell_median_comp_cov', 'cell_mean_comp_cov', 'cell_std_comp_cov',
        'cell_nonzero_comp_cov', 'cell_max_streak_comp_cov', 'cell_min_streak_comp_cov',
        'cell_min_junc_cov', 'cell_max_junc_cov', 'cell_median_junc_cov', 'cell_mean_junc_cov', 'cell_std_junc_cov',
        'cell_nonzero_junc_cov', 'cell_max_streak_junc_cov', 'cell_min_streak_junc_cov',
        'label_general', 'label_specific'
    ]
    df = pd.DataFrame(all_data, columns=feature_names)
    
    # Merge with meta_prediction
    meta_columns = [
        'data_source', 'sample_size', 'transcript_id',
        'bottleneck_coverage', 'highest_coverage', 
        'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov',
        'y_prob_general'
    ]

    #meta_columns = ['data_source', 'sample_size', 'transcript_id', 'bottleneck_coverage', 'highest_coverage', 'y_prob_general']
    df = pd.merge(df, meta_prediction[meta_columns], 
                  on=['data_source', 'sample_size', 'transcript_id'], 
                  how='left')
    
    # Statistics of cell-specific df
    general_pos_total = np.sum(df['label_general'])
    general_neg_total = len(df) - general_pos_total
    specific_pos_total = np.sum(df['label_specific'])
    specific_neg_total = len(df) - specific_pos_total

    # Count positives and negatives for y_prob_general < 0.2
    low_prob_df = df[df['y_prob_general'] < 0.2]
    general_pos = np.sum(low_prob_df['label_general'])
    general_neg = len(low_prob_df) - general_pos
    specific_pos = np.sum(low_prob_df['label_specific'])
    specific_neg = len(low_prob_df) - specific_pos
    
    print(f"For y_prob_general < 0.2:")
    print(f"label_general: {general_pos}/{general_pos_total} positives, {general_neg}/{general_neg_total} negatives")
    print(f"label_specific: {specific_pos}/{specific_pos_total} positives, {specific_neg}/{specific_neg_total} negatives")
    
    df = df[df['y_prob_general'] >= 0.2]

    # Split data into train and test sets
    train_df = df[df['chr'].isin(train_chrm)]
    test_df = df[~df['chr'].isin(train_chrm)]
    
    return train_df, test_df

def process_and_save_specific_data_noSpecificGT(data_dir, data_source, sample_size_list, train_chrm, meta_prediction):
    """Process and save data for all sample sizes."""
    for sample_size in sample_size_list:
        # Filter meta_prediction for current sample_size
        meta_prediction_filtered = meta_prediction[meta_prediction['sample_size'] == sample_size]
        
        train_df, test_df = combine_general_to_specific_noSpecificGT(data_dir, data_source, sample_size, train_chrm, meta_prediction_filtered)
        
        os.makedirs(f'{TRAIN_DIR}/train_test_data', exist_ok=True)
        
        train_df.to_csv(f'{TRAIN_DIR}/train_test_data/{data_source}_{sample_size}_train.csv', index=False)
        test_df.to_csv(f'{TRAIN_DIR}/train_test_data/{data_source}_{sample_size}_test.csv', index=False)
        
        print(f"Saved {data_source} data for sample size {sample_size}")

# %%
def load_and_combine_specific_data(train_dir, data_sources):
    """Load, print statistics, and combine train and test data from CSV files."""
    train_data = []
    test_data = []
    # Define column types, specifying 'chr' as string
    dtype_dict = {'chr': str}
    
    print("\nData Statistics:")
    for source in data_sources:
        if source == 'HEK293T':
            sample_sizes = HUMAN_SAMPLE_SIZE_LIST
        elif source == 'Mouse-Fibroblast':
            sample_sizes = MOUSE_SAMPLE_SIZE_LIST
        else:
            raise ValueError(f"Unknown data source: {source}")
        
        for size in sample_sizes:
            train_file = f'{train_dir}/train_test_data/{source}_{size}_train.csv'
            test_file = f'{train_dir}/train_test_data/{source}_{size}_test.csv'
            
            if not os.path.exists(train_file) or not os.path.exists(test_file):
                print(f"Warning: Files not found for {source} - Sample size {size}")
                continue
            
            train_df = pd.read_csv(train_file, dtype=dtype_dict)
            test_df = pd.read_csv(test_file, dtype=dtype_dict)
            #train_df = train_df[train_df['y_prob_general']>=0.2]
            #test_df = test_df[test_df['y_prob_general']>=0.2]
            train_data.append(train_df)
            test_data.append(test_df)
            
            print(f"\n{source} - Sample size {size}:")
            
            # Calculate statistics for train data
            train_total = len(train_df)
            train_general_pos = train_df['label_general'].sum()
            train_specific_pos = train_df['label_specific'].sum()
            
            # Calculate statistics for test data
            test_total = len(test_df)
            test_general_pos = test_df['label_general'].sum()
            test_specific_pos = test_df['label_specific'].sum()
            
            # Print statistics
            print(f"Train - Total: {train_total}, General Positive: {train_general_pos}, Specific Positive: {train_specific_pos}")
            print(f"Test - Total: {test_total}, General Positive: {test_general_pos}, Specific Positive: {test_specific_pos}")
            
            # Calculate and print precisions
            total_instances = train_total + test_total
            total_general_pos = train_general_pos + test_general_pos
            total_specific_pos = train_specific_pos + test_specific_pos
            
            print(f"General Precision: {total_general_pos}/{total_instances} = {total_general_pos/total_instances:.3f}")
            print(f"Specific Precision: {total_specific_pos}/{total_instances} = {total_specific_pos/total_instances:.3f}")
    
    combined_train = pd.concat(train_data, ignore_index=True)
    combined_test = pd.concat(test_data, ignore_index=True)
    
    print(f"\nTotal train data shape: {combined_train.shape}")
    print(f"Total test data shape: {combined_test.shape}")
    
    # Print overall statistics
    total_instances = len(combined_train) + len(combined_test)
    total_general_pos = combined_train['label_general'].sum() + combined_test['label_general'].sum()
    total_specific_pos = combined_train['label_specific'].sum() + combined_test['label_specific'].sum()
    
    print(f"\nOverall Statistics:")
    print(f"Total Instances: {total_instances}")
    print(f"Total General Positives: {total_general_pos}")
    print(f"Total Specific Positives: {total_specific_pos}")
    print(f"Overall General Precision: {total_general_pos}/{total_instances} = {total_general_pos/total_instances:.3f}")
    print(f"Overall Specific Precision: {total_specific_pos}/{total_instances} = {total_specific_pos/total_instances:.3f}")
    
    return combined_train, combined_test

# %%
def process_and_analyze_specific_data(data_dir, data_source, sample_size_list, train_chrm, meta_prediction):
    """Process, analyze, and return combined data for all sample sizes."""
    train_data = []
    test_data = []
    
    print(f"\nData Statistics for {data_source}:")
    
    for sample_size in sample_size_list:
        # Filter meta_prediction for current sample_size
        meta_prediction_filtered = meta_prediction[meta_prediction['sample_size'] == sample_size]
        
        # Process data
        train_df, test_df = combine_general_to_specific(data_dir, data_source, sample_size, train_chrm, meta_prediction_filtered)
        
        # train_df = train_df[train_df['y_prob_general'] >= 0.2]
        # test_df = test_df[test_df['y_prob_general'] >= 0.2]
        
        train_data.append(train_df)
        test_data.append(test_df)
        
        # Calculate and print statistics
        print(f"Sample size {sample_size}:")
        print_meta_specific_label_statistics(train_df, test_df)
    
    # Combine all data
    combined_train = pd.concat(train_data, ignore_index=True)
    combined_test = pd.concat(test_data, ignore_index=True)
    
    # Print overall statistics
    print("\nOverall Statistics:")
    print_overall_statistics(combined_train, combined_test)
    
    return combined_train, combined_test

def print_meta_specific_label_statistics(train_df, test_df):
    """Print statistics for train and test data."""

    # Calculate statistics for train data
    train_total = len(train_df)
    train_general_pos = train_df['label_general'].sum()
    train_specific_pos = train_df['label_specific'].sum()
    
    # Calculate statistics for test data
    test_total = len(test_df)
    test_general_pos = test_df['label_general'].sum()
    test_specific_pos = test_df['label_specific'].sum()
    
    print(f"Train - Total: {train_total}, General Positive: {train_general_pos}, Specific Positive: {train_specific_pos}")
    print(f"Test - Total: {test_total}, General Positive: {test_general_pos}, Specific Positive: {test_specific_pos}")
    
    total_instances = train_total + test_total
    total_general_pos = train_general_pos + test_general_pos
    total_specific_pos = train_specific_pos + test_specific_pos
    
    print(f"General Precision: {total_general_pos}/{total_instances} = {total_general_pos/total_instances:.3f}")
    print(f"Specific Precision: {total_specific_pos}/{total_instances} = {total_specific_pos/total_instances:.3f}\n")

def print_overall_statistics(combined_train, combined_test):
    """Print overall statistics for combined data."""
    total_instances = len(combined_train) + len(combined_test)
    total_general_pos = combined_train['label_general'].sum() + combined_test['label_general'].sum()
    total_specific_pos = combined_train['label_specific'].sum() + combined_test['label_specific'].sum()
    
    print(f"Total Instances: {total_instances}")
    print(f"Total General Positives: {total_general_pos}")
    print(f"Total Specific Positives: {total_specific_pos}")
    print(f"Overall General Precision: {total_general_pos}/{total_instances} = {total_general_pos/total_instances:.3f}")
    print(f"Overall Specific Precision: {total_specific_pos}/{total_instances} = {total_specific_pos/total_instances:.3f}")

# %%
def train_and_evaluate_specific_model(X_train, y_train, X_test, y_test, model_name):
    """Train a Random Forest model, evaluate its performance, and attach probabilities to test data."""
    feature_columns = [
        'sample_size',
        'input_genes', 'input_transcripts', 'output_transcripts',
        'cell_supporting_junctions', 'cell_ratio_supporting_junctions',
        'cell_min_comp_cov', 'cell_max_comp_cov', 'cell_median_comp_cov', 'cell_mean_comp_cov', 'cell_std_comp_cov',
        'cell_nonzero_comp_cov', 'cell_max_streak_comp_cov', 'cell_min_streak_comp_cov',
        'cell_min_junc_cov', 'cell_max_junc_cov', 'cell_median_junc_cov', 'cell_mean_junc_cov', 'cell_std_junc_cov',
        'cell_nonzero_junc_cov', 'cell_max_streak_junc_cov', 'cell_min_streak_junc_cov', 
        'bottleneck_coverage', 'highest_coverage', 
        'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov',
        #'y_prob_general'
    ]
    
    X_train = X_train[feature_columns]
    X_test = X_test[feature_columns]

    class_weights = {0: 1, 1: 5}
    model = RandomForestClassifier(n_estimators=100, max_depth=12, random_state=42, class_weight=class_weights)
    model.fit(X_train, y_train)
    
    # Train performance
    y_train_pred = model.predict(X_train)
    y_train_prob = model.predict_proba(X_train)[:, 1]
    train_precision = precision_score(y_train, y_train_pred)
    train_recall = recall_score(y_train, y_train_pred)
    train_roc_auc = roc_auc_score(y_train, y_train_prob)
    
    print(f'{model_name} Train Results:')
    print(f'Train Precision: {train_precision:.3f}')
    print(f'Train Recall: {train_recall:.3f}')
    print(f'Train ROC curve (area = {train_roc_auc:.3f})')
    
    # Test performance
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    test_precision = precision_score(y_test, y_pred)
    test_recall = recall_score(y_test, y_pred)
    test_roc_auc = roc_auc_score(y_test, y_prob)
    
    print(f'\n{model_name} Test Results:')
    print(f'Test Precision: {test_precision:.3f}')
    print(f'Test Recall: {test_recall:.3f}')
    print(f'Test ROC curve (area = {test_roc_auc:.3f})')
    
    # Save the model
    model_filename = f'../models/{model_name}_roc={test_roc_auc:.3f}.joblib'
    joblib.dump(model, model_filename)
    
    return model, model_filename, y_prob

def load_and_evaluate_specific_model(model_filename, X_test, y_test):
    """Load a saved model, evaluate its performance on test data, and attach probabilities."""
    model = joblib.load(model_filename)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    test_precision = precision_score(y_test, y_pred)
    test_recall = recall_score(y_test, y_pred)
    test_roc_auc = roc_auc_score(y_test, y_prob)
    
    print(f'Loaded Model Results:')
    print(f'Test Precision: {test_precision:.3f}')
    print(f'Test Recall: {test_recall:.3f}')
    print(f'Test ROC AUC: {test_roc_auc:.3f}')
    
    return y_prob

# %%
def process_individual_cells(combined_test, data_dir, data_source, sample_size):
    """
    Process individual cell data using precomputed predictions and save results.
    """
    cell_dir = f"{data_dir}/merge_{sample_size}/linkmerge_{sample_size}_sgtf"
    output_dir = f"{TRAIN_DIR}/predictions/{data_source}/{sample_size}"
    os.makedirs(output_dir, exist_ok=True)
    
    for cell_id in range(1, sample_size + 1):
        # Filter combined_test for the current data_source and sample_size
        current_test = combined_test[(combined_test['data_source'] == data_source) & 
                                     (combined_test['sample_size'] == sample_size) &
                                     (combined_test['cell_id'] == cell_id)]
        
        tmap_file = os.path.join(cell_dir, f"individual.{cell_id}.{cell_id}.gtf.tmap")
        output_file = os.path.join(output_dir, f"{data_source}.{sample_size}.cell-{cell_id}.csv")
        
        if not os.path.exists(tmap_file):
            print(f"Warning: tmap file for cell {cell_id} in {data_source}, sample size {sample_size} does not exist.")
            continue
        
        # Read tmap file
        tmap_data = []
        with open(tmap_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                fields = line.strip().split('\t')
                ref_gene_id, ref_id, class_code, qry_id = fields[0], fields[1], fields[2], fields[4]
                label = 1 if class_code == '=' else 0
                tmap_data.append([ref_gene_id, ref_id, class_code, qry_id, label])
        
        tmap_df = pd.DataFrame(tmap_data, columns=['ref_gene_id', 'ref_id', 'class_code', 'qry_id', 'label'])
        
        # Get intersection of transcript_ids
        cell_transcripts = set(tmap_df['qry_id'])
        test_transcripts = set(current_test['transcript_id'])
        common_transcripts = cell_transcripts.intersection(test_transcripts)
        
        # Filter tmap_df and current_test for common transcripts
        tmap_df = tmap_df[tmap_df['qry_id'].isin(common_transcripts)]
        cell_test = current_test[current_test['transcript_id'].isin(common_transcripts)]
        
        # Merge predictions with tmap data
        merged_df = pd.merge(tmap_df, cell_test[['transcript_id', 'y_prob_general', 'y_prob', 'label_specific']], 
                             left_on='qry_id', right_on='transcript_id')
        
        # Sort by y_prob in descending order
        merged_df = merged_df.sort_values('y_prob', ascending=False).reset_index(drop=True)
        
        # Calculate cumulative matches and precision
        merged_df['match_cum'] = merged_df['label_specific'].cumsum()
        merged_df['precision_cum'] = merged_df['match_cum'] / (merged_df.index + 1)
        
        # Select and order columns
        result_df = merged_df[['ref_gene_id', 'ref_id', 'class_code', 'qry_id', 'y_prob_general', 'label', 'label_specific', 'y_prob', 'match_cum', 'precision_cum']]
        result_df = result_df.rename(columns={'label': 'label_general'})
        
        # Save results
        result_df.to_csv(output_file, index=False)
        
        # Print statistics
        total_instances = len(result_df)
        positive_cases = result_df['label_specific'].sum()
        negative_cases = total_instances - positive_cases
        print(f"Cell {cell_id} in {data_source}, sample size {sample_size}:")
        #print(f"  Total instances: {total_instances}")
        #print(f"  Positive cases: {positive_cases}")
        #print(f"  Negative cases: {negative_cases}")

# %%
def compare_with_other_tools(data_source, sample_size_list, tools, predictions_dir, plot_data_dir):
    """
    Compare beaver predictions with other tools and generate a comparison table.
    """
    comparison_results = []
    plot_data = {tool: {'mat': [], 'pre': [], 'beaver_mat': [], 'beaver_pre': []} for tool in tools}
    os.makedirs(plot_data_dir, exist_ok=True)

    for sample_size in sample_size_list:
        sample_result = {'Data Source': data_source, 'Sample Size': sample_size}
        
        for tool in tools:
            # Read tool data from CSV
            tool_csv_path = os.path.join(TRAIN_BASE_DIR, "realDataset-tools-otherchrs", tool, f"{data_source}_statistics.csv")
            tool_data = pd.read_csv(tool_csv_path)
            tool_data_sample = tool_data[tool_data['sample_size'] == sample_size]
            
            tool_mat = tool_data_sample['matching_specific'].tolist()
            tool_pre = tool_data_sample['precision_specific'].tolist()
            
            beaver_adjusted_pre = []
            beaver_adjusted_mat = []
            beaver_actual_mat = []
            
            csv_data = []

            for cell_id in range(1, sample_size + 1):
                prediction_file = os.path.join(predictions_dir, str(sample_size), f"{data_source}.{sample_size}.cell-{cell_id}.csv")
                
                if not os.path.exists(prediction_file):
                    print(f"Warning: Prediction file for cell {cell_id} in {data_source}, sample size {sample_size} does not exist.")
                    continue
                
                cell_predictions = pd.read_csv(prediction_file)
                
                adjusted_mat = tool_mat[cell_id - 1]
                tool_precision = tool_pre[cell_id - 1]
                matching_indices = np.where(cell_predictions['match_cum'].values == adjusted_mat)[0]
                
                if len(matching_indices) > 0:
                    adjusted_index = matching_indices[0]
                else:
                    # Find the closest precision
                    precision_diff = cell_predictions['precision_cum'] - tool_precision
                    lower_precision_indices = np.where(precision_diff <= 0)[0]
                    
                    if len(lower_precision_indices) > 0:
                        adjusted_index = lower_precision_indices[0]
                    else:
                        # If no lower precision found, use the last row
                        adjusted_index = len(cell_predictions) - 1
                    
                    print(f"Warning: No exact match found for tool {tool} in cell {cell_id}. Using closest precision.")

                adjusted_precision = cell_predictions.loc[adjusted_index, 'precision_cum']
                beaver_mat = cell_predictions.loc[adjusted_index, 'match_cum']
                
                beaver_adjusted_pre.append(adjusted_precision)
                beaver_adjusted_mat.append(adjusted_mat)
                beaver_actual_mat.append(beaver_mat)
                
                csv_data.append({
                    'cell_id': cell_id,
                    'adjusted_mat': adjusted_mat,
                    'beaver_mat': beaver_mat,
                    f'{tool}_precision': round(tool_precision*100, 1),
                    'beaver_precision': round(adjusted_precision*100, 1)
                })
            
            if beaver_adjusted_pre:
                beaver_avg_precision = np.mean(beaver_adjusted_pre)* 100
                tool_avg_precision = np.mean(tool_pre) * 100
                delta = (beaver_avg_precision - tool_avg_precision) #/ tool_avg_precision * 100
                
                sample_result[f'{tool} #matching'] = round(np.mean(tool_mat), 1)
                sample_result[f'{tool} average precision'] = round(tool_avg_precision, 1)
                sample_result[f'beaver vs {tool} precision'] = round(beaver_avg_precision, 1)
                sample_result[f'{tool} delta'] = round(delta, 1)
                
                plot_data[tool]['mat'].append(tool_mat)
                plot_data[tool]['pre'].append(tool_pre)
                plot_data[tool]['beaver_mat'].append(beaver_actual_mat)
                plot_data[tool]['beaver_pre'].append(beaver_adjusted_pre)
                
                # Save CSV for this tool and sample size
                csv_filename = os.path.join(plot_data_dir, f"{data_source}_{tool}_{sample_size}.csv")
                pd.DataFrame(csv_data).to_csv(csv_filename, index=False)
                #print(f"Saved plot data to {csv_filename}")
        
        comparison_results.append(sample_result)
        #print(f"Processed sample size {sample_size} for {data_source}")
        print(sample_result)
    
    return pd.DataFrame(comparison_results), plot_data

# %%
def plot_scatter_separate(sample_size_list, figDir, plot_data_dir, data_source):
    tool_colors = {
        'beaver': '#FF5733',
        'transmeta': '#3498DB',
        'psiclass': '#8B668B',
        'aletsch': '#C87C78',
        'stringtie2': '#F39C12',
        'scallop2': '#397125'
    }

    tool_names = {
        'beaver': 'Beaver',
        'transmeta': 'TransMeta',
        'psiclass': 'PsiCLASS',
        'aletsch': 'Aletsch',
        'stringtie2': 'StringTie2',
        'scallop2': 'Scallop2'
    }

    os.makedirs(figDir, exist_ok=True)

    n_samples = len(sample_size_list)
    n_cols = 3
    n_rows = math.ceil(n_samples / n_cols)

    for tool in tool_colors.keys():
        if tool == 'beaver':
            continue  # Skip beaver

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5*n_rows), squeeze=False)
        #fig.suptitle(f"Comparison of {tool} and beaver - {data_source}", fontsize=24)
        
        for i, sample_size in enumerate(sample_size_list):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col]
            
            csv_filename = os.path.join(plot_data_dir, f"{data_source}_{tool}_{sample_size}.csv")
            if not os.path.exists(csv_filename):
                print(f"Warning: CSV file {csv_filename} does not exist.")
                continue
            
            df = pd.read_csv(csv_filename)
            
            ax.scatter(df[f'{tool}_precision'], df['adjusted_mat'], 
                       color=tool_colors[tool], label=tool_names[tool], 
                       alpha=0.6, s=40, marker='o')
            ax.scatter(df['beaver_precision'], df['beaver_mat'], 
                       color=tool_colors['beaver'], label='Beaver', 
                       alpha=0.6, s=40, marker='o')

            ax.set_xlabel("Precision (%)", fontsize=20)
            ax.set_ylabel("# Matching Transcripts", fontsize=20)
            ax.set_title(f"# Cells = {sample_size}", fontsize=20)
            
            ax.tick_params(axis='both', which='major', labelsize=14)

            x_min = min(df[f'{tool}_precision'].min(), df['beaver_precision'].min())
            x_max = max(df[f'{tool}_precision'].max(), df['beaver_precision'].max())
            y_min = df['adjusted_mat'].min()
            y_max = df['adjusted_mat'].max()
            
            x_padding = (x_max - x_min) * 0.05
            y_padding = (y_max - y_min) * 0.05
            
            ax.set_xlim(max(0, x_min - x_padding), min(100, x_max + x_padding))
            ax.set_ylim(max(0, y_min - y_padding), y_max + y_padding)

            ax.grid(True, axis='y', linestyle='--', alpha=0.7)
            ax.set_axisbelow(True)

        for i in range(n_samples, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            fig.delaxes(axes[row, col])

        plt.tight_layout(rect=[0, 0.06, 1, 0.98])  # Adjusted to leave more space for legend
        
        # Create a single legend for the entire figure
        handles, labels = ax.get_legend_handles_labels()
        leg = fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.02), 
                 ncol=2, fontsize=22, markerscale=3, handlelength=10,
                 columnspacing=10, labelspacing=0.3, handletextpad=0.1, frameon=False)
        
        figname = os.path.join(figDir, f'{data_source}_{tool}_comparison.pdf')
        plt.savefig(figname, format="pdf", bbox_inches="tight")
        plt.show()
        plt.close()

    print(f"Plots saved in: {figDir}")

# %%
def process_meta_models(data_source, load_size_list, sample_size_list):
    print(f"Processing Meta Models for {data_source}:")
    
    beaver_dir = HUMAN_BEAVER_DIR if data_source == 'HEK293T' else MOUSE_BEAVER_DIR
    
    process_and_save_meta_data(beaver_dir, data_source, load_size_list, TRAIN_CHRM)
    meta_train, meta_test = load_and_combine_meta_data(TRAIN_DIR, [data_source])
    meta_model, meta_model_filename, meta_yprob, meta_train_yprob = train_and_evaluate_meta_model(
        meta_train.drop('label_general', axis=1),
        meta_train['label_general'],
        meta_test.drop('label_general', axis=1),
        meta_test['label_general'],
        f'model_real{data_source}_{BEAVER}_meta'
    )
    print(f"{data_source} Meta Model saved at: {meta_model_filename}")
    meta_train['y_prob_general'] = meta_train_yprob
    meta_test['y_prob_general'] = meta_yprob

    # Process meta predictions for individual cells
    print(f"\nSave meta(general) predictions for individual cells ({data_source}):")
    for sample_size in sample_size_list:
        process_meta_individual_cells(meta_test, beaver_dir, data_source, sample_size)

    return meta_train, meta_test

def process_mixed_meta_model(data_sources):
    print("\nProcessing Mixed Meta Model:")
    
    meta_mix_train, meta_mix_test = load_and_combine_meta_data(TRAIN_DIR, data_sources)
    mix_model, mix_model_filename, mix_yprob, mix_train_yprob = train_and_evaluate_meta_model(
        meta_mix_train.drop('label_general', axis=1),
        meta_mix_train['label_general'],
        meta_mix_test.drop('label_general', axis=1),
        meta_mix_test['label_general'],
        f'model_realMix_meta_{BEAVER}'
    )
    print(f"Mixed Model saved at: {mix_model_filename}")
    meta_mix_train['y_prob_general'] = mix_train_yprob
    meta_mix_test['y_prob_general'] = mix_yprob

    return meta_mix_train, meta_mix_test

def process_cell_specific_models(data_source, load_size_list, sample_size_list, meta_train, meta_test):
    print(f"\nProcessing Cell-specific Models for {data_source}:")
    
    beaver_dir = HUMAN_BEAVER_DIR if data_source == 'HEK293T' else MOUSE_BEAVER_DIR
    
    combined_meta = pd.concat([meta_train, meta_test], ignore_index=True)
    process_and_save_specific_data_noSpecificGT(beaver_dir, data_source, load_size_list, TRAIN_CHRM, combined_meta)
    combined_train, combined_test = load_and_combine_specific_data(TRAIN_DIR, [data_source])

    #combined_train, combined_test = process_and_analyze_specific_data(
        #beaver_dir, data_source, sample_size_list, TRAIN_CHRM, combined_meta
        #)

    model, model_filename, yprob = train_and_evaluate_specific_model(
        combined_train.drop('label_specific', axis=1),
        combined_train['label_specific'],
        combined_test.drop('label_specific', axis=1),
        combined_test['label_specific'],
        f'model_real{data_source}_{BEAVER}_specific'
    )
    print(f"{data_source} Cell-specific Model saved at: {model_filename}")
    combined_test['y_prob'] = yprob

    # Process individual cells
    print(f"\nProcessing individual cells for {data_source}:")
    for sample_size in sample_size_list:
        process_individual_cells(combined_test, beaver_dir, data_source, sample_size)

def evaluate_and_save_comparison(data_source, sample_size_list):
    print(f"\nEvaluating and Plotting Results for {data_source}:")
    
    comparison, plot_data = compare_with_other_tools(data_source, 
                                                     sample_size_list,
                                                     ['transmeta', 'psiclass', 'aletsch', 'stringtie2', 'scallop2'], 
                                                     os.path.join(TRAIN_DIR, "predictions", data_source), 
                                                     os.path.join(TRAIN_DIR, "plots"))
    comparison.to_csv(os.path.join(TRAIN_DIR, f"{data_source}_comparison_results.csv"), index=False)
    
# %%
#process_all_tools('HEK293T')
#process_all_tools('Mouse-Fibroblast')
#process_all_tools_other_chrs('HEK293T')
#process_all_tools_other_chrs('Mouse-Fibroblast')

# %%
HUMAN_LOAD_SIZE_LIST = [2, 5, 10, 20, 30, 50, 80, 100, 150, 192]
#HUMAN_LOAD_SIZE_LIST = [2]
MOUSE_LOAD_SIZE_LIST = [2, 5, 10, 30, 50, 80, 100, 150, 200, 300, 369]
#MOUSE_LOAD_SIZE_LIST = [2]

# %%
# 1. Meta Model Stage
meta_human_train, meta_human_test = process_meta_models('HEK293T', HUMAN_LOAD_SIZE_LIST, HUMAN_SAMPLE_SIZE_LIST)
meta_mouse_train, meta_mouse_test = process_meta_models('Mouse-Fibroblast', MOUSE_LOAD_SIZE_LIST, MOUSE_SAMPLE_SIZE_LIST)

#meta_mix_train, meta_mix_test = process_mixed_meta_model(['HEK293T', 'Mouse-Fibroblast'])
#meta_human_train = meta_mix_train[meta_mix_train['data_source'] == 'HEK293T']
#meta_human_test = meta_mix_test[meta_mix_test['data_source'] == 'HEK293T']
#meta_mouse_train = meta_mix_train[meta_mix_train['data_source'] == 'Mouse-Fibroblast']
#meta_mouse_test = meta_mix_test[meta_mix_test['data_source'] == 'Mouse-Fibroblast']

# %%
# 2. Cell-specific Model Stage
process_cell_specific_models('HEK293T', HUMAN_LOAD_SIZE_LIST, HUMAN_SAMPLE_SIZE_LIST, meta_human_train, meta_human_test)
process_cell_specific_models('Mouse-Fibroblast', MOUSE_LOAD_SIZE_LIST, MOUSE_SAMPLE_SIZE_LIST, meta_mouse_train, meta_mouse_test)

# %%
# 3. Evaluation and Plotting Stage
evaluate_and_save_comparison('HEK293T', HUMAN_SAMPLE_SIZE_LIST)
evaluate_and_save_comparison('Mouse-Fibroblast', MOUSE_SAMPLE_SIZE_LIST)

# %%
human_fig_dir = os.path.join(TRAIN_DIR, "plots", 'HEK293T')
os.makedirs(human_fig_dir, exist_ok=True)
HUMAN_SAMPLE_SIZE_PLOT_LIST = [5, 10, 30, 50, 100, 192]
print("\nGenerating plots for human data:")
plot_scatter_separate(HUMAN_SAMPLE_SIZE_PLOT_LIST, human_fig_dir, os.path.join(TRAIN_DIR, "plots"), 'HEK293T')

mouse_fig_dir = os.path.join(TRAIN_DIR, "plots", 'Mouse-Fibroblast')
os.makedirs(mouse_fig_dir, exist_ok=True)
MOUSE_SAMPLE_SIZE_PLOT_LIST = [10, 30, 50, 100, 200, 369]
print("\nGenerating plots for mouse data:")
plot_scatter_separate(MOUSE_SAMPLE_SIZE_PLOT_LIST, mouse_fig_dir, os.path.join(TRAIN_DIR, "plots"), 'Mouse-Fibroblast')
    

# %%
