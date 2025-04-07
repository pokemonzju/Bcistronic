###ORF CSV File Format Requirements
This script requires an input CSV file (orf_csv_file) with the following mandatory columns:
Trans_ID: Transcript ID.
ORF_ID: ORF information.

ORF_ID Format
The ORF_ID values must follow this specific format:
ORF<Number>_<Transcript_ID>:<Start>:<End>[;...]
Multiple ORFs can be separated by semicolons ;.

Example
Trans_ID        ORF_ID
FBti0018861     ORF8_FBti0018861:5590:6960;ORF9_FBti0018861:2507:5575  
FBti0018862     ORF10_FBti0018862:2813:5575;ORF6_FBti0018862:5590:6960  
FBti0018875     ORF1_FBti0018875:2244:4964;ORF2_FBti0018875:340:2247
