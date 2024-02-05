import os 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import matplotlib.patheffects as path_effects

def add_median_labels(ax, fmt='.3f'):
    lines = ax.get_lines()
    boxes = [c for c in ax.get_children() if type(c).__name__ == 'PathPatch']
    lines_per_box = int(len(lines) / len(boxes))
    for median in lines[4:len(lines):lines_per_box]:
        x, y = (data.mean() for data in median.get_data())
        # choose value depending on horizontal or vertical plot orientation
        value = x if (median.get_xdata()[1] - median.get_xdata()[0]) == 0 else y
        text = ax.text(x, y, f'{value:{fmt}}', ha='center', va='center',
                       fontweight='bold', color='white')
        # create median-colored border around white text for contrast
        text.set_path_effects([
            path_effects.Stroke(linewidth=3, foreground=median.get_color()),
            path_effects.Normal(),
        ])

def plot_param(gait_data_p1, trial, comp_trial, param):

    ya_data = gait_data_p1[(gait_data_p1['group'] == 'YA') & (gait_data_p1['trial_type'] == trial)][param].dropna()
    oa_data = gait_data_p1[(gait_data_p1['group'] == 'OA') & (gait_data_p1['trial_type'] == trial)][param].dropna()
    pd_data = gait_data_p1[(gait_data_p1['group'] == 'PD') & (gait_data_p1['trial_type'] == trial)][param].dropna()

    print(ya_data.mean())
    print(oa_data.mean())
    print(pd_data.mean())
    print(oa_data.mean() - ya_data.mean())
    print(pd_data.mean() - oa_data.mean())

    print(f"\nAll blocks: {param}")
    t, p = scipy.stats.ttest_ind(oa_data[0:len(ya_data)], ya_data)
    print('OA vs. YA t = ', str(t), ' p = ', str(p))
    t, p = scipy.stats.ttest_ind(pd_data, oa_data[0:len(pd_data)])
    print('PD vs. OA t = ', str(t), ' p = ', str(p))

    ya_data_dt = gait_data_p1[(gait_data_p1['group'] == 'YA') & (gait_data_p1['trial_type'] == comp_trial)][param].dropna()
    oa_data_dt = gait_data_p1[(gait_data_p1['group'] == 'OA') & (gait_data_p1['trial_type'] == comp_trial)][param].dropna()
    pd_data_dt = gait_data_p1[(gait_data_p1['group'] == 'PD') & (gait_data_p1['trial_type'] == comp_trial)][param].dropna()
    
    t, p = scipy.stats.ttest_ind(ya_data_dt, ya_data)
    print('YA DT vs. YA ST t = ', str(t), ' p = ', str(p))
    t, p = scipy.stats.ttest_ind(oa_data_dt, oa_data)
    print('OA DT vs. OA ST t = ', str(t), ' p = ', str(p))
    t, p = scipy.stats.ttest_ind(pd_data_dt, pd_data)
    print('PD DT vs. PD ST t = ', str(t), ' p = ', str(p))

    box_plot = sns.boxplot(data=gait_data_p1, x=param, y="group", hue="trial_type")
    add_median_labels(box_plot)
    box_plot.set_axisbelow(True)
    box_plot.xaxis.grid(True) 
    plt.show()

if __name__ == "__main__":

    pd.set_option('display.max_rows', None)

    param_file = "../Data/imu_gait_parameters.csv"
    gait_param_data = pd.read_csv(param_file)

    ya_ids = pd.read_csv("../Data/identifiers_YA.csv")['id_nummer'].to_list()
    oa_ids = pd.read_csv("../Data/identifiers_OA.csv")['id_nummer'].to_list()
    pd_ids = pd.read_csv("../Data/identifiers_PD.csv")['id_nummer'].to_list()

    gait_data_p1 = gait_param_data[(gait_param_data['session'] == 'protocol1')]
    gait_data_p1['Walking Speed LR'] = gait_data_p1[["Walking Speed R", "Walking Speed L"]].mean(axis=1)
    gait_data_p1['Step Time LR'] = gait_data_p1[["Step Time R", "Step Time L"]].mean(axis=1)

    ids = sorted(list(set(gait_data_p1['subject'].to_list())))
    for id in ya_ids:
        gait_data_p1.loc[gait_data_p1.subject == id, 'group'] = "YA"
    for id in oa_ids:
        gait_data_p1.loc[gait_data_p1.subject == id, 'group'] = "OA"
    for id in pd_ids:
        gait_data_p1.loc[gait_data_p1.subject == id, 'group'] = "PD"

    param_file = "../Data/imu_variability_parameters.csv"
    variability_param_data = pd.read_csv(param_file)
    variability_param_data_p1 = variability_param_data[(variability_param_data['session'] == 'protocol1')]

    # Plot comparisons    
    trial = "Straight_walking"
    comp_trial = "Straight_walking_and_Aud_Stroop"

    # Get condition variables
    avg_frames = []
    for id in ids:
        st_data = gait_data_p1[(gait_data_p1['subject'] == id) & (gait_data_p1['trial_type'] == trial)]
        dt_data = gait_data_p1[(gait_data_p1['subject'] == id) & (gait_data_p1['trial_type'] == comp_trial)]
        avg_ws_st = st_data["Walking Speed LR"].mean(axis=0)
        avg_ws_dt = dt_data["Walking Speed LR"].mean(axis=0)
        dt_diff = avg_ws_dt - avg_ws_st
        dt_cost_walk_speed = -(dt_diff / avg_ws_st) * 100 # DTE from Kelly et al 2010

        var_data = variability_param_data_p1[(variability_param_data_p1['subject'] == id) & 
                                                     (variability_param_data_p1['trial_type'] == trial)]
        st_step_time_var = var_data['Step Time Variability']
        if st_step_time_var.empty:
            st_step_time_var = np.NaN
        else:
            st_step_time_var = st_step_time_var.values[0]

        avg_data = {
            'subject': [id],
            'group': [st_data['group'].values[0]],
            'dt_cost_walk_speed': [dt_cost_walk_speed],
            'st_walk_speed': [avg_ws_st],
            'st_step_time_var': [st_step_time_var]
        }
        avg_frame = pd.DataFrame(data=avg_data)
        avg_frames.append(avg_frame)
    avg_frame = pd.concat(avg_frames)

    print(avg_frame)

    # Plot 
    box_plot = sns.boxplot(data=avg_frame, x="dt_cost_walk_speed", y="group")
    add_median_labels(box_plot)
    box_plot.set_axisbelow(True)
    box_plot.xaxis.grid(True) 
    plt.show()

    box_plot = sns.boxplot(data=avg_frame, x="st_walk_speed", y="group")
    add_median_labels(box_plot)
    box_plot.set_axisbelow(True)
    box_plot.xaxis.grid(True) 
    plt.show()

    box_plot = sns.boxplot(data=avg_frame, x="st_step_time_var", y="group")
    add_median_labels(box_plot)
    box_plot.set_axisbelow(True)
    box_plot.xaxis.grid(True) 
    plt.show()

    avg_frame.to_csv("../Data/validation_data_imu_gait_parameters.csv", index=False)
    