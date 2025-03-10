using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;
using System.Text;
using System.Numerics;
using System.Collections.Generic;
using UnityEditor;

public class strucmap_create_usecsv : MonoBehaviour
{
    static public int size=64;

    void Start()
    {
        Texture3D texture;
        readfile();
        //Debug.Log(Refl1(30,1,1.5,1.6,200,20));
        texture=create_3dtex();
        AssetDatabase.CreateAsset(texture, "Assets/LUTs/strucmap_m3d8-N8-new2.asset");//変える
    }

    // Update is called once per frame
    void Update()
    {
        
    }
    public double[,,,] strucTable = new double[size,size,size,3];//縦、横、色(*3しているのはRGBあるため)
    void readfile()
    {
        string filePath_1=@"C:\Users\長谷川　紘大\Documents\長谷川\Unity\Squid_HDRP\Assets\csv\Multilayer_m3d8-N8-new2.csv";//変える
        
        // StreamReaderクラスをインスタンス化
        StreamReader reader = new StreamReader(filePath_1, Encoding.GetEncoding("UTF-8"));
        
        // 最後まで読み込む
        int count=0;
        int zdep=0;
        while (reader.Peek() >= 0)
        {
            // 読み込んだ文字列をカンマ区切りで配列に格納
            string[] cols = reader.ReadLine().Split(',');
              // 表示
            //Debug.Log(cols.Length);
            string vol;
            for(int i=0;i<size*3;i++){
                vol=cols[i];
                if(vol=="nan"){
                    vol="0";
                }
                float a=float.Parse(vol);
                strucTable[i/3,count,zdep,i%3]=a/255;//正規化をする
                //Debug.Log(a);
            }
            
            count++;
            if(count==size){
                count=0;
                zdep++;
            }
        }
        reader.Close();
    }
    
    Texture3D create_3dtex(){
        Texture3D texture;
        Color[] colorArray = new Color[size * size * size];
        texture = new Texture3D (size, size, size, TextureFormat.RGBA32, true);
        float r = 1.0f / (size-1.0f);
        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++) {
                for (int z = 0; z < size; z++) {
                    //Color c = new Color (x*r, y*r, z*r, 1.0f);
                    //Color c = new Color (0.0f,1.0f,0.0f, 1.0f);
                    float R=(float)strucTable[x,y,z,0];//青
                    float G=(float)strucTable[x,y,z,1];//緑
                    float B=(float)strucTable[x,y,z,2];//赤
                    Color c = new Color (R, G, B, 1.0f);
                    //colorArray[y + (x * size) + (z * size * size)] = c;
                    colorArray[x + (y * size) + (z * size * size)] = c;
                    //Debug.Log(x*r);
                }
            }
        }
        texture.SetPixels (colorArray);
        texture.Apply ();
        return texture;
    }
    
}
