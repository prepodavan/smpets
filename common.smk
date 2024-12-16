def sample_kind_from_id(wildcards):
    for sample in config["samples"]["dna"]["ids"]:
        if sample in wildcards.sample:
            return 'dna'
    for sample in config["samples"]["rna"]["ids"]:
        if sample in wildcards.sample:
            return 'rna'
    return None

def calling_filters(wildcards):
    return config[sample_kind_from_id(wildcards)]['filters']['calling']

def trimmomatic_se_filters(wildcards):
    return config[sample_kind_from_id(wildcards)]['filters']['trimmomatic_se']

def trimmomatic_pe_filters(wildcards):
    return config[sample_kind_from_id(wildcards)]['filters']['trimmomatic_pe']

