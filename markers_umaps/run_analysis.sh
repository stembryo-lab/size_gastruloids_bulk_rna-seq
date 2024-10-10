#!/bin/bash

dataset2analyze=$1

if [[ $dataset2analyze == 'prisca' ]]
then
    SAMPLES=(
    S1-0h
    S2-48h
    S3-52h
    S4-56h
    S6-96h
    S7.1-108h
    S7.2-108h
    S8.2-120h
    SBR-24h
    SBR-36h
    SBR-48h
    SBR-60h
    SBR-72h
    SBR-84h_A
    SBR-84h_B)
fi

if [[ $dataset2analyze == 'pijuan' ]]
then
    SAMPLES=(
    pijuan_E6.5_18
    pijuan_E6.5_1
    pijuan_E6.5_5
    pijuan_E6.75_7
    pijuan_E7.0_10
    pijuan_E7.0_14
    pijuan_E7.0_15
    pijuan_E7.0_30
    pijuan_E7.0_31
    pijuan_E7.0_32
    pijuan_E7.25_23
    pijuan_E7.25_26
    pijuan_E7.25_27
    pijuan_E7.5_19
    pijuan_E7.5_20
    pijuan_E7.5_2
    pijuan_E7.5_3
    pijuan_E7.5_4
    pijuan_E7.5_6
    pijuan_E7.75_12
    pijuan_E7.75_13
    pijuan_E7.75_8
    pijuan_E7.75_9
    pijuan_E8.0_16
    pijuan_E8.0_33
    pijuan_E8.0_34
    pijuan_E8.0_35
    pijuan_E8.25_24
    pijuan_E8.25_25
    pijuan_E8.25_28
    pijuan_E8.5_17
    pijuan_E8.5_29
    pijuan_E8.5_36
    pijuan_E8.5_37
    # pijuan_E8.5_86
    # pijuan_E8.5_87
    # pijuan_E8.5_88
    # pijuan_E8.5_89
    # pijuan_E8.5_90
    # pijuan_E8.5_91
    # pijuan_E8.5_92
    # pijuan_E8.5_93
    # pijuan_E8.5_94
    # pijuan_E8.5_95
    # pijuan_E8.5_96
    # pijuan_E8.5_97
    pijuan_E8.75_54
    pijuan_E8.75_55
    pijuan_E8.75_56
    pijuan_E8.75_57
    pijuan_E8.75_58
    pijuan_E8.75_59
    # pijuan_E8.75_60
    # pijuan_E8.75_61
    # pijuan_E8.75_62
    # pijuan_E8.75_63
    # pijuan_E8.75_64
    # pijuan_E8.75_65
    # pijuan_E8.75_66
    # pijuan_E8.75_67
    # pijuan_E8.75_68
    # pijuan_E8.75_69
    pijuan_E9.0_38
    pijuan_E9.0_39
    pijuan_E9.0_40
    pijuan_E9.0_41
    pijuan_E9.0_42
    pijuan_E9.0_43
    # pijuan_E9.0_44
    # pijuan_E9.0_45
    # pijuan_E9.0_46
    # pijuan_E9.0_47
    # pijuan_E9.0_48
    # pijuan_E9.0_49
    # pijuan_E9.0_50
    # pijuan_E9.0_51
    # pijuan_E9.0_52
    # pijuan_E9.0_53
    pijuan_E9.25_70
    pijuan_E9.25_71
    pijuan_E9.25_72
    pijuan_E9.25_73
    pijuan_E9.25_74
    pijuan_E9.25_75
    # pijuan_E9.25_76
    # pijuan_E9.25_77
    # pijuan_E9.25_78
    # pijuan_E9.25_79
    # pijuan_E9.25_80
    # pijuan_E9.25_81
    # pijuan_E9.25_82
    # pijuan_E9.25_83
    # pijuan_E9.25_84
    # pijuan_E9.25_85
    pijuan_E9.5_100
    pijuan_E9.5_101
    pijuan_E9.5_102
    pijuan_E9.5_103
    pijuan_E9.5_104
    pijuan_E9.5_105
    # pijuan_E9.5_98
    # pijuan_E9.5_99
    )
fi



# RUN QC AND ANALYSIS
for sample in ${SAMPLES[@]}
do
    echo "Start analyzing ${sample} ..."
    python analysis_by_batch.py $sample data/mouse tiff
    echo "${sample} analysis completed"
    
done

echo "All samples processed"


# RUN INTEGRATION
echo 'Integrate prisca samples'

#python analysis_integration.py

echo 'Analysis fininshed!'



