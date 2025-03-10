using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ScoreManager : MonoBehaviour {

    public GameObject score_object = null; // Textオブジェクト
    public float score_num = 0f; // スコア変数

      // 初期化
      void Start () {
      }

      // 更新
      void Update () {
        // オブジェクトからTextコンポーネントを取得
        Text score_text = score_object.GetComponent<Text> ();

        score_num = Mathf.Cos((Time.time) * Mathf.PI / 7);
        score_num = 0.5f * score_num + 0.5f;


        // テキストの表示を入れ替える
        score_text.text = "D=" + score_num.ToString("F2");
      }
}