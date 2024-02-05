import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == "__main__":

    #pd.set_option('display.max_rows', None)

    demo_file = "../Data/basic_demographics.csv"
    demo_data = pd.read_csv(demo_file)

    crf_file = "../Data/REDcap_data/CRF_data.csv"
    crf_data = pd.read_csv(crf_file)

    updrs_file = "../Data/REDcap_data/UPDRS_data.csv"
    updrs_data = pd.read_csv(updrs_file)

    ya_ids = pd.read_csv("../Data/identifiers_YA.csv")['id_nummer'].to_list()
    oa_ids = pd.read_csv("../Data/identifiers_OA.csv")['id_nummer'].to_list()
    pd_ids = pd.read_csv("../Data/identifiers_PD.csv")['id_nummer'].to_list()

    ids = sorted(list(set(demo_data['subject'].to_list())))
    for id in ya_ids:
        demo_data.loc[demo_data.subject == id, 'group'] = "YA"
    for id in oa_ids:
        demo_data.loc[demo_data.subject == id, 'group'] = "OA"
    for id in pd_ids:
        demo_data.loc[demo_data.subject == id, 'group'] = "PD"
    print(demo_data)

    groups = ['YA', 'OA', 'PD']
    i = 0
    for group_ids in [ya_ids, oa_ids, pd_ids]:

        print('Group: %s' % (groups[i]))

        demo_data_group = demo_data[demo_data['subject'].isin(group_ids)]
        crf_data_group = crf_data[crf_data['id_nummer'].isin(group_ids)]

        age_group = demo_data_group['age']
        sex_group = demo_data_group['sex']
        edu_group = crf_data_group['crf_utbildning_ar']
        fd_group = demo_data_group['frandin_grimby']
        print("Group age - Range %d-%d / mean %f / SD %f " % (age_group.min(), age_group.max(), age_group.mean(), age_group.std())) 
        print("Group edu - Range %d-%d / mean %f / SD %f " % (edu_group.min(), edu_group.max(), edu_group.mean(), edu_group.std())) 
        print("Group sex - Female n %d / percent %f" % ((sex_group == 1).sum(), (sex_group == 1).sum() / len(sex_group))) 
        print("Group frd - Range %d-%d / mean %f / SD %f " % (fd_group.min(), fd_group.max(), fd_group.mean(), fd_group.std())) 

        i +=1

    # Barplot for HY
    updrs_data_pd = updrs_data[updrs_data['id_nummer'].isin(pd_ids)]
    ax = sns.countplot(data=updrs_data_pd, x="mdsupdrs3_hy")
    ax.bar_label(ax.containers[0])
    plt.show()