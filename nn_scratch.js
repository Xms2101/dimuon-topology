    //(0,0), (1,1)
    //(X,1), (X,X)
    //(Y,1), (Y,Y)
    //(Z,1), (Z,Z)
    //(1,0), (2,0), (X,0), (Y,0), (Z,0)
    
    
    
    if ((neutronsA_count == 0 && neutronsC_count == 0)){ //(0,0)
        h_a_n[0]     ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[0]   ->Fill(kperp);
        h_k_lin[0]   ->Fill(kperp);
        h_mass[0]    ->Fill(mass);
        h_rapidity[0]->Fill(rapidity);
        h_acop[0]    ->Fill(acop);
        if (abs(acop) < 0.008){
            h_acop_s[0]->Fill(acop);
            h_a_n_s[0] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[0] ->Fill(kperp);
            h_mass_s[0]    ->Fill(mass);
            h_rapidity_s[0]->Fill(rapidity);
        }
    }
    if ((neutronsA_count == 1 && neutronsC_count == 1)){ //(1,1)
        h_a_n[1]     ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[1]   ->Fill(kperp);
        h_k_lin[1]   ->Fill(kperp);
        h_mass[1]    ->Fill(mass);
        h_rapidity[1]->Fill(rapidity);
        h_acop[1]    ->Fill(acop);
        if (abs(acop) < 0.008){
            h_acop_s[1]->Fill(acop);
            h_a_n_s[1] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[1] ->Fill(kperp);
            h_mass_s[1]    ->Fill(mass);
            h_rapidity_s[1]->Fill(rapidity);
        }
    }    
    if ((neutronsA_count > 0 && neutronsC_count == 1) || (neutronsA_count == 1 && neutronsC_count > 0)){ //(X,1)
        h_a_n[2]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[2]   ->Fill(kperp);
        h_k_lin[2]   ->Fill(kperp);
        h_mass[2]    ->Fill(mass);
        h_rapidity[2]->Fill(rapidity);
        //h_acop_n3      ->Fill(acop);
        h_acop[2]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n3->Fill(acop);
            h_acop_s[2]->Fill(acop);
            h_a_n_s[2] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[2] ->Fill(kperp);
            h_mass_s[2]    ->Fill(mass);
            h_rapidity_s[2]->Fill(rapidity);
        }
    }
    if (neutronsA_count > 0 && neutronsC_count > 0){ //(X,X)
        h_a_n[3]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[3]   ->Fill(kperp);
        h_k_lin[3]   ->Fill(kperp);
        h_mass[3]    ->Fill(mass);
        h_rapidity[3]->Fill(rapidity);
        //h_acop_n4      ->Fill(acop);
        h_acop[3]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n4->Fill(acop);
            h_acop_s[3]->Fill(acop);
            h_a_n_s[3] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[3] ->Fill(kperp);
            h_mass_s[3]    ->Fill(mass);
            h_rapidity_s[3]->Fill(rapidity);
        }
    }
    if ((neutronsA_count > 1 && neutronsC_count == 1) || (neutronsA_count == 1 && neutronsC_count > 1)  && (neutronsA_count <= 5) && (neutronsC_count <= 5)){ //(Y,1)
        h_a_n[4]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[4]   ->Fill(kperp);
        h_k_lin[4]   ->Fill(kperp);
        h_mass[4]    ->Fill(mass);
        h_rapidity[4]->Fill(rapidity);
        //h_acop_n5      ->Fill(acop);
        h_acop[4]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n5->Fill(acop);
            h_acop_s[4]->Fill(acop);
            h_a_n_s[4] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[4] ->Fill(kperp);
            h_mass_s[4]    ->Fill(mass);
            h_rapidity_s[4]->Fill(rapidity);
        }
    }
    if (neutronsA_count > 1 && neutronsC_count > 1  && (neutronsA_count <= 5) && (neutronsC_count <= 5)){ //(Y,Y)
        h_a_n[5]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[5]   ->Fill(kperp);
        h_k_lin[5]   ->Fill(kperp);
        h_mass[5]    ->Fill(mass);
        h_rapidity[5]->Fill(rapidity);
        //h_acop_n6      ->Fill(acop);
        h_acop[5]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n6->Fill(acop);
            h_acop_s[5]->Fill(acop);
            h_a_n_s[5] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[5] ->Fill(kperp);
            h_mass_s[5]    ->Fill(mass);
            h_rapidity_s[5]->Fill(rapidity);
        }
    }
    if ((neutronsA_count > 5 && neutronsC_count == 1) || (neutronsA_count == 1 && neutronsC_count > 5)){ //(Z,1)
        h_a_n[6]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[6]   ->Fill(kperp);
        h_k_lin[6]   ->Fill(kperp);
        h_mass[6]    ->Fill(mass);
        h_rapidity[6]->Fill(rapidity);
        //h_acop_n7      ->Fill(acop);
        h_acop[6]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n7->Fill(acop);
            h_acop_s[6]->Fill(acop);
            h_a_n_s[6] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[6] ->Fill(kperp);
            h_mass_s[6]    ->Fill(mass);
            h_rapidity_s[6]->Fill(rapidity);
        }
    }
    if (neutronsA_count > 5 && neutronsC_count > 5){ //(Z,Z)
        h_a_n[7]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[7]   ->Fill(kperp);
        h_k_lin[7]   ->Fill(kperp);
        h_mass[7]    ->Fill(mass);
        h_rapidity[7]->Fill(rapidity);
        //h_acop_n8      ->Fill(acop);
        h_acop[7]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n8->Fill(acop);
            h_acop_s[7]->Fill(acop);
            h_a_n_s[7] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[7] ->Fill(kperp);
            h_mass_s[7]    ->Fill(mass);
            h_rapidity_s[7]->Fill(rapidity);
        }    
    }
    if ((neutronsA_count == 1 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count == 1)){ //(1,0)
        h_a_n[8]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[8]   ->Fill(kperp);
        h_k_lin[8]   ->Fill(kperp);
        h_mass[8]    ->Fill(mass);
        h_rapidity[8]->Fill(rapidity);
        //h_acop_n9      ->Fill(acop);
        h_acop[8]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n9->Fill(acop);
            h_acop_s[8]->Fill(acop);
            h_a_n_s[8] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[8] ->Fill(kperp);
            h_mass_s[8]    ->Fill(mass);
            h_rapidity_s[8]->Fill(rapidity);
        }  
    }
    if ((neutronsA_count == 2 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count == 2)){ //(2,0)
        h_a_n[9]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[9]   ->Fill(kperp);
        h_k_lin[9]   ->Fill(kperp);
        h_mass[9]    ->Fill(mass);
        h_rapidity[9]->Fill(rapidity);
        //h_acop_n10      ->Fill(acop);
        h_acop[9]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n10->Fill(acop);
            h_acop_s[9]->Fill(acop);
            h_a_n_s[9] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[9] ->Fill(kperp);
            h_mass_s[9]    ->Fill(mass);
            h_rapidity_s[9]->Fill(rapidity);
        } 
    }
    if ((neutronsA_count == 0 && neutronsC_count > 0) || (neutronsA_count > 0 && neutronsC_count == 0)){ //(X,0)
        h_a_n[10]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[10]   ->Fill(kperp);
        h_k_lin[10]   ->Fill(kperp);
        h_mass[10]    ->Fill(mass);
        h_rapidity[10]->Fill(rapidity);
        //h_acop_n11      ->Fill(acop);
        h_acop[10]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n11->Fill(acop);
            h_acop_s[10]->Fill(acop);
            h_a_n_s[10] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[10] ->Fill(kperp);
            h_mass_s[10]    ->Fill(mass);
            h_rapidity_s[10]->Fill(rapidity);
        }  
    }
    if ((neutronsA_count > 1 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count > 1)  && (neutronsA_count <= 5) && (neutronsC_count <= 5)){ //(Y,0)
        h_a_n[11]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[11]   ->Fill(kperp);
        h_k_lin[11]   ->Fill(kperp);
        h_mass[11]    ->Fill(mass);
        h_rapidity[11]->Fill(rapidity);
        //h_acop_n12      ->Fill(acop);
        h_acop[11]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n12->Fill(acop);
            h_acop_s[11]->Fill(acop);
            h_a_n_s[11] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[11] ->Fill(kperp);
            h_mass_s[11]    ->Fill(mass);
            h_rapidity_s[11]->Fill(rapidity);
        }
    }
    if ((neutronsA_count > 5 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count > 5)){ //(Z,0)
        h_a_n[12]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[12]   ->Fill(kperp);
        h_k_lin[12]   ->Fill(kperp);
        h_mass[12]    ->Fill(mass);
        h_rapidity[12]->Fill(rapidity);
        //h_acop_n10      ->Fill(acop);
        h_acop[12]    ->Fill(acop);
        if (abs(acop) < 0.008){
            //h_acop_s_n10->Fill(acop);
            h_acop_s[12]->Fill(acop);
            h_a_n_s[12] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[12] ->Fill(kperp);
            h_mass_s[12]    ->Fill(mass);
            h_rapidity_s[12]->Fill(rapidity);
        } 
    }