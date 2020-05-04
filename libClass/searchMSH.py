# ==========================================
# Code created by Leandro Marques at 02/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code searches a file until 12 directories


# ------------------------------------------------------------------
# Use:
# directory = search_file.Find(mesh_name)
# if directory == 'File does not found':
#  sys.exit()
# ------------------------------------------------------------------




import os

def Find(_file):
  file = _file
  
  # Exit any directory to /home/user
  path = os.getcwd()
  os.chdir(path)

  # Directories list only
  files1 = filter(os.path.isdir,os.listdir(path))

  breaking = 0

  if file in os.listdir(path):
   path = path

  # If there is any directories, access them
  else:   
   for dir1 in files1:
    
    # Access directory 1
    path = path + '/' + dir1
    os.chdir(path)
    
    # Directories list only
    files2 = filter(os.path.isdir,os.listdir(path))
    
    if file in os.listdir(path):
     breaking = 1
     break

    # There are not directories
    elif files2 == []:
     aa = path.split('/')
     directory = ''
   
     # Reduce a directory
     for i in range(1,len(aa)-1):
      bb = aa[i]
      directory = directory + '/' + bb
 
     path = directory
     os.chdir(path)
   
    # If there is any directories, access them
    # The statements onwards are equivalent to previous statement
    else:
     for dir2 in files2:
      path = path + '/' + dir2
      os.chdir(path)
      files3 = filter(os.path.isdir,os.listdir(path))
      
      if file in os.listdir(path):
       breaking = 1
       break
 
      elif files3 == []:
       aa = path.split('/')
       directory = ''
    
       for i in range(1,len(aa)-1):
        bb = aa[i]
        directory = directory + '/' + bb
 
       path = directory
       os.chdir(path)
   
      else:
       for dir3 in files3:
        path = path + '/' + dir3
        os.chdir(path)
        files4 = filter(os.path.isdir,os.listdir(path))
        
        if file in os.listdir(path):
         breaking = 1
         break
 
        elif files4 == []:
         aa = path.split('/')
         directory = ''
    
         for i in range(1,len(aa)-1):
          bb = aa[i]
          directory = directory + '/' + bb
 
         path = directory
         os.chdir(path)
  
        else:
         for dir4 in files4:
          path = path + '/' + dir4
          os.chdir(path)
          files5 = filter(os.path.isdir,os.listdir(path))
          
          if file in os.listdir(path):
           breaking = 1
           break
 
          elif files5 == []:
           aa = path.split('/')
           directory = ''
    
           for i in range(1,len(aa)-1):
            bb = aa[i]
            directory = directory + '/' + bb
 
           path = directory
           os.chdir(path)
   
          else:
           for dir5 in files5:
            path = path + '/' + dir5
            os.chdir(path)
            files6 = filter(os.path.isdir,os.listdir(path))
            
            if file in os.listdir(path):
             breaking = 1
             break
 
            elif files6 == []:
             aa = path.split('/')
             directory = ''
    
             for i in range(1,len(aa)-1):
              bb = aa[i]
              directory = directory + '/' + bb
 
             path = directory
             os.chdir(path)
  
            else:
             for dir6 in files6:
              path = path + '/' + dir6
              os.chdir(path)
              files7 = filter(os.path.isdir,os.listdir(path))
              
              if file in os.listdir(path):
               breaking = 1
               break
 
              elif files7 == []:
               aa = path.split('/')
               directory = ''
    
               for i in range(1,len(aa)-1):
                bb = aa[i]
                directory = directory + '/' + bb
 
               path = directory
               os.chdir(path)
  
              else:
               for dir7 in files7:
                path = path + '/' + dir7
                os.chdir(path)
                files8 = filter(os.path.isdir,os.listdir(path))
                
                if file in os.listdir(path):
                 breaking = 1
                 break
 
                elif files8 == []:
                 aa = path.split('/')
                 directory = ''
    
                 for i in range(1,len(aa)-1):
                  bb = aa[i]
                  directory = directory + '/' + bb
 
                 path = directory
                 os.chdir(path)
  
                else:
                 for dir8 in files8:
                  path = path + '/' + dir8
                  os.chdir(path)
                  files9 = filter(os.path.isdir,os.listdir(path))
                  
                  if file in os.listdir(path):
                   breaking = 1
                   break
 
                  elif files9 == []:
                   aa = path.split('/')
                   directory = ''
    
                   for i in range(1,len(aa)-1):
                    bb = aa[i]
                    directory = directory + '/' + bb
 
                   path = directory
                   os.chdir(path)
  

                  else:
                   for dir9 in files9:
                    path = path + '/' + dir9
                    os.chdir(path)
                    files10 = filter(os.path.isdir,os.listdir(path))
                    
                    if file in os.listdir(path):
                     breaking = 1
                     break
  
                    elif files10 == []:
                     aa = path.split('/')
                     directory = ''
    
                     for i in range(1,len(aa)-1):
                      bb = aa[i]
                      directory = directory + '/' + bb
  
                     path = directory
                     os.chdir(path)
   
                    else:
                     for dir10 in files10:
                      path = path + '/' + dir10
                      os.chdir(path)
                      files11 = filter(os.path.isdir,os.listdir(path))
                      
                      if file in os.listdir(path):
                       breaking = 1
                       break
  
                      elif files11 == []:
                       aa = path.split('/')
                       directory = ''
     
                       for i in range(1,len(aa)-1):
                        bb = aa[i]
                        directory = directory + '/' + bb
  
                       path = directory
                       os.chdir(path)
   
                      else:
                       for dir11 in files11:
                        path = path + '/' + dir11
                        os.chdir(path)
                        files12 = filter(os.path.isdir,os.listdir(path))
                        
                        if file in os.listdir(path):
                         breaking = 1
                         break
 
                        elif files12 == []:
                         aa = path.split('/')
                         directory = ''
    
                         for i in range(1,len(aa)-1):
                          bb = aa[i]
                          directory = directory + '/' + bb
  
                         path = directory
                         os.chdir(path)
  

                        else:
                         aa = path.split('/')
                         directory = ''
     
                         for i in range(1,len(aa)-1):
                          bb = aa[i]
                          directory = directory + '/' + bb
  
                         path = directory
                         os.chdir(path)
  
  
                       if breaking == 1:
                        break
                       
                       # The search final 
                       else:  
                        aa = path.split('/')
                        directory = ''
    
                        for i in range(1,len(aa)-1):
                         bb = aa[i]
                         directory = directory + '/' + bb
 
                        path = directory
                        os.chdir(path)
 
                     # Break directories for-loop 
                     if breaking == 1:
                       break
                     
                     # If file does not found, reduce a directory
                     # The statements onwards are equivalent to previous statement
                     else:  
                      aa = path.split('/')
                      directory = ''
    
                      for i in range(1,len(aa)-1):
                       bb = aa[i]
                       directory = directory + '/' + bb
 
                      path = directory
                      os.chdir(path)
 
                   if breaking == 1:
                    break
    
                   else:  
                    aa = path.split('/')
                    directory = ''
    
                    for i in range(1,len(aa)-1):
                     bb = aa[i]
                     directory = directory + '/' + bb
 
                    path = directory
                    os.chdir(path)
 
                 if breaking == 1:
                  break
    
                 else:  
                  aa = path.split('/')
                  directory = ''
    
                  for i in range(1,len(aa)-1):
                   bb = aa[i]
                   directory = directory + '/' + bb
 
                  path = directory
                  os.chdir(path)
 
               if breaking == 1:
                break
    
               else:  
                aa = path.split('/')
                directory = ''
    
                for i in range(1,len(aa)-1):
                 bb = aa[i]
                 directory = directory + '/' + bb
 
                path = directory
                os.chdir(path)
 
             if breaking == 1:
              break
   
             else:  
              aa = path.split('/')
              directory = ''
    
              for i in range(1,len(aa)-1):
               bb = aa[i]
               directory = directory + '/' + bb

              path = directory
              os.chdir(path)
 
           if breaking == 1:
            break
   
           else:  
            aa = path.split('/')
            directory = ''
   
            for i in range(1,len(aa)-1):
             bb = aa[i]
             directory = directory + '/' + bb

            path = directory
            os.chdir(path)

         if breaking == 1:
          break
   
         else:  
          aa = path.split('/')
          directory = ''
   
          for i in range(1,len(aa)-1):
           bb = aa[i]
           directory = directory + '/' + bb

          path = directory
          os.chdir(path)

       if breaking == 1:
        break
   
       else:  
        aa = path.split('/')
        directory = ''
   
        for i in range(1,len(aa)-1):
         bb = aa[i]
         directory = directory + '/' + bb

        path = directory
        os.chdir(path)

     if breaking == 1:
      break
   
     else:  
      aa = path.split('/')
      directory = ''
   
      for i in range(1,len(aa)-1):
       bb = aa[i]
       directory = directory + '/' + bb

      path = directory
      os.chdir(path)

  if breaking == 1:
   return path
   
  else:
   print ""
   print " Error: File not found"
   print ""
   
   path = 'File not found'
   return path

