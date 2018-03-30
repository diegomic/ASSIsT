def snp_classes():
    #return ['Recovered_2Het','2ClusterNull', '4ClusterNull', 'OneClassMissing',
    #        'OneHomoRare_HWE','OneHomoRare_NotHWE','Robust', 'ShiftedHomo',
    #        'NullAllele', 'DistortedAndUnexpSegreg', 'FalseSNP','Failed', 'NonInformative']
    return ['AB_2_sub-clusters',           # 'Recovered_2Het'
            'Null_2_Cluster',  # '2ClusterNull'
            'Null_4_Cluster',  # '4ClusterNull'
            'OneClassMissing',            # 'OneClassMissing'
            'OneHomozygRare_HWE',      # 'OneHomoRare_HWE'
            'OneHomozygRare_NotHWE',  # 'OneHomoRare_NotHWE'
            'Robust',                       # 'Robust'
            'ShiftedHomo',                  # 'ShiftedHomo'
            'NullAllele-Failed',            # 'NullAllele'
            'DistortedAndUnexpSegreg',                    # 'DistortedAndUnexpSegreg'
            'Monomorphic',                    # 'FalseSNP'
            'Failed',                       # 'Failed'
            'Not-Informative']#             # 'NonInformative'

def classes_approved_F1F2_snps():
    #===============================================================================================
    # exported_class = ['2ClusterNull','4ClusterNull','Recovered_2Het','Robust','OneHomoRare_HWE']
    #===============================================================================================       
    return ['Robust',                       # 'Robust'
            'AB_2_sub-clusters',           # 'Recovered_2Het'
            'Null_2_Cluster',  # '2ClusterNull'
            'Null_4_Cluster']  # '4ClusterNull'
            #'OneHomozygRare_HWE',      # 'OneHomoRare_HWE'
            #'OneHomozygRare_NotHWE']  # 'OneHomoRare_NotHWE'
            

def classes_approved_Germplasm_snps():
    #===============================================================================================
    # exported_class = ['DistortedAndUnexpSegreg', 'OneHomoRare_HWE', 
    #                           'OneHomoRare_NotHWE','Robust']
    #===============================================================================================       
    return ['Robust',                       # 'Robust'
            'OneHomozygRare_HWE',      # 'OneHomoRare_HWE'
            'OneHomozygRare_NotHWE',  # 'OneHomoRare_NotHWE'
            'DistortedAndUnexpSegreg']                    # 'DistortedAndUnexpSegreg'

def classes_discarded_F1F2_snps():
    #===============================================================================================
    # non_exported_class = ['Failed','NullAllele', 'FalseSNP','NonInformative', 'OneClassMissing', 'ShiftedHomo']
    #===============================================================================================
    return ['Monomorphic',                    # 'FalseSNP'
            'Failed',                       # 'Failed'
            'DistortedAndUnexpSegreg',                    # 'DistortedAndUnexpSegreg'
            'NullAllele-Failed',            # 'NullAllele'
            'OneClassMissing',              # 'ShiftedHomo'
            'NotInformative']              # 'NonInformative'


def classes_discarded_Germplasm_snps():
    #===============================================================================================
    # non_exported_class = ['Failed','NullAllele', 'FalseSNP','NonInformative', 'OneClassMissing', 'ShiftedHomo']
    #===============================================================================================
    return ['Monomorphic',                    # 'FalseSNP'
            'Failed',                       # 'Failed'
            'ShiftedHomo',                  # 'ShiftedHomo'
            'NullAllele-Failed',            # 'NullAllele'
            'NotInformative']#             # 'NonInformative'